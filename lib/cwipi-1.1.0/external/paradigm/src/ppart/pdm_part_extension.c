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
#include "pdm_distant_neighbor.h"
#include "pdm_logging.h"
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_order.h"
#include "pdm_error.h"
#include "pdm_part_extension.h"
#include "pdm_part_to_part.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_part_extension_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_vtk.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

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

static
PDM_g_num_t*
_compute_offset_ln_to_gn_by_domain
(
  int              n_domain,
  int             *n_part,
  int            **pn_entity,
  PDM_g_num_t   ***pentity_ln_to_gn,
  PDM_MPI_Comm     comm
)
{

  PDM_g_num_t *shift_by_domain_loc = PDM_array_const_gnum(n_domain, 0);
  PDM_g_num_t *shift_by_domain     = (PDM_g_num_t *) malloc((n_domain+1) * sizeof(PDM_g_num_t));

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {

      int          _pn_entity        = pn_entity       [i_domain][i_part];
      PDM_g_num_t *_pentity_ln_to_gn = pentity_ln_to_gn[i_domain][i_part];
      for(int i = 0; i < _pn_entity; ++i) {
        shift_by_domain_loc[i_domain] = PDM_MAX(shift_by_domain_loc[i_domain], _pentity_ln_to_gn[i]);
      }
    }
  }

  shift_by_domain[0] = 0;
  PDM_MPI_Allreduce(shift_by_domain_loc, &shift_by_domain[1], n_domain, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_array_accumulate_gnum(shift_by_domain, n_domain+1);

  free(shift_by_domain_loc);

  return shift_by_domain;
}


static
void
_offset_ln_to_gn_by_domain
(
  int              n_domain,
  int             *n_part,
  int            **pn_entity,
  PDM_g_num_t   ***pentity_ln_to_gn,
  PDM_g_num_t     *shift_by_domain,
  int              sens
)
{
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      int          _pn_entity        = pn_entity       [i_domain][i_part];
      PDM_g_num_t *_pentity_ln_to_gn = pentity_ln_to_gn[i_domain][i_part];
      for(int i = 0; i < _pn_entity; ++i) {
        _pentity_ln_to_gn[i] = _pentity_ln_to_gn[i] + sens * shift_by_domain[i_domain];
      }
    }
  }
}

static
void
_offset_parts_by_domain
(
  PDM_part_extension_t *part_ext,
  int                   sens
)
{
  int **pn_cell       = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_face       = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_edge       = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_vtx        = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_face_group = malloc(part_ext->n_domain * sizeof(int *));

  PDM_g_num_t ***cell_ln_to_gn       = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***face_ln_to_gn       = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***edge_ln_to_gn       = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***vtx_ln_to_gn        = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***face_group_ln_to_gn = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    pn_cell      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_face      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_edge      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_vtx       [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_face_group[i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));

    cell_ln_to_gn      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    face_ln_to_gn      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    edge_ln_to_gn      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    vtx_ln_to_gn       [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    face_group_ln_to_gn[i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      pn_cell            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_cell;
      pn_face            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
      pn_edge            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
      pn_vtx             [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_face_group      [i_domain][i_part] = 0;
      if (part_ext->parts[i_domain][i_part].n_face_group > 0) {
        pn_face_group      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_bound_idx[part_ext->parts[i_domain][i_part].n_face_group];
      }

      cell_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].cell_ln_to_gn;
      face_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
      edge_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
      vtx_ln_to_gn       [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      face_group_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_group_ln_to_gn;

    }
  }

  // Here we go
  if(sens == 1) {
    assert(part_ext->shift_by_domain_cell == NULL);
    part_ext->shift_by_domain_cell = _compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                        part_ext->n_part,
                                                                        pn_cell,
                                                                        cell_ln_to_gn,
                                                                        part_ext->comm);

    assert(part_ext->shift_by_domain_face == NULL);
    part_ext->shift_by_domain_face = _compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                        part_ext->n_part,
                                                                        pn_face,
                                                                        face_ln_to_gn,
                                                                        part_ext->comm);

    assert(part_ext->shift_by_domain_edge == NULL);
    part_ext->shift_by_domain_edge = _compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                        part_ext->n_part,
                                                                        pn_edge,
                                                                        edge_ln_to_gn,
                                                                        part_ext->comm);

    assert(part_ext->shift_by_domain_vtx == NULL);
    part_ext->shift_by_domain_vtx = _compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                       part_ext->n_part,
                                                                       pn_vtx,
                                                                       vtx_ln_to_gn,
                                                                       part_ext->comm);

    assert(part_ext->shift_by_domain_face_group == NULL);
    part_ext->shift_by_domain_face_group = _compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                       part_ext->n_part,
                                                                       pn_face_group,
                                                                       face_group_ln_to_gn,
                                                                       part_ext->comm);
  }


  _offset_ln_to_gn_by_domain(part_ext->n_domain,
                             part_ext->n_part,
                             pn_cell,
                             cell_ln_to_gn,
                             part_ext->shift_by_domain_cell,
                             sens);


  _offset_ln_to_gn_by_domain(part_ext->n_domain,
                             part_ext->n_part,
                             pn_face,
                             face_ln_to_gn,
                             part_ext->shift_by_domain_face,
                             sens);

  _offset_ln_to_gn_by_domain(part_ext->n_domain,
                             part_ext->n_part,
                             pn_edge,
                             edge_ln_to_gn,
                             part_ext->shift_by_domain_edge,
                             sens);

  _offset_ln_to_gn_by_domain(part_ext->n_domain,
                             part_ext->n_part,
                             pn_vtx,
                             vtx_ln_to_gn,
                             part_ext->shift_by_domain_vtx,
                             sens);

  _offset_ln_to_gn_by_domain(part_ext->n_domain,
                             part_ext->n_part,
                             pn_face_group,
                             face_group_ln_to_gn,
                             part_ext->shift_by_domain_face_group,
                             sens);


  // Attention il faut shifter tout les border_ln_to_gn aussi et deduire le border_i_domain
  // (Recherche dicotomique dans le shift by domain)
  if(part_ext->pdi !=NULL) {
    part_ext->n_interface = PDM_part_domain_interface_n_interface_get(part_ext->pdi);
  }

  if(part_ext->pdi != NULL && part_ext->n_domain > 1 && part_ext->n_interface > 0) {
    abort();
  }

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(pn_cell      [i_domain]);
    free(pn_face      [i_domain]);
    free(pn_edge      [i_domain]);
    free(pn_vtx       [i_domain]);
    free(pn_face_group[i_domain]);

    free(cell_ln_to_gn      [i_domain]);
    free(face_ln_to_gn      [i_domain]);
    free(edge_ln_to_gn      [i_domain]);
    free(vtx_ln_to_gn       [i_domain]);
    free(face_group_ln_to_gn[i_domain]);
  }

  free(pn_cell      );
  free(pn_face      );
  free(pn_edge      );
  free(pn_vtx       );
  free(pn_face_group);

  free(cell_ln_to_gn      );
  free(face_ln_to_gn      );
  free(edge_ln_to_gn      );
  free(vtx_ln_to_gn       );
  free(face_group_ln_to_gn);
}

static
void
_shift_ln_to_gn
(
  int          n_entity,
  PDM_g_num_t *entity_ln_to_gn,
  PDM_g_num_t  shift,
  int          sens
)
{
  for(int i = 0; i < n_entity; ++i) {
    entity_ln_to_gn[i] = entity_ln_to_gn[i] + sens * shift;
  }
}


static
void
_offset_results_by_domain
(
  PDM_part_extension_t *part_ext
)
{
  PDM_UNUSED(part_ext);
  // You need to shift border_* and fill the array border_dom_id

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell        = part_ext->parts[i_domain][i_part].n_cell;
      int n_border_cell = part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][n_cell];
      PDM_g_num_t *cell_ln_to_gn = part_ext->border_cell_ln_to_gn[shift_part+i_part];
      _shift_ln_to_gn(n_border_cell, cell_ln_to_gn, part_ext->shift_by_domain_cell[i_domain], -1);

      int n_face        = part_ext->parts[i_domain][i_part].n_face;
      int n_border_face = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      PDM_g_num_t *face_ln_to_gn = part_ext->border_face_ln_to_gn[shift_part+i_part];
      _shift_ln_to_gn(n_border_face, face_ln_to_gn, part_ext->shift_by_domain_face[i_domain], -1);

      int n_edge        = part_ext->parts[i_domain][i_part].n_edge;
      if(n_edge > 0 && part_ext->edge_edge_extended_idx != NULL ) {
        int n_border_edge = part_ext->edge_edge_extended_idx[shift_part+i_part][n_edge];
        PDM_g_num_t *edge_ln_to_gn = part_ext->border_edge_ln_to_gn[shift_part+i_part];
        _shift_ln_to_gn(n_border_edge, edge_ln_to_gn, part_ext->shift_by_domain_edge[i_domain], -1);
      }

      int n_vtx        = part_ext->parts[i_domain][i_part].n_vtx;
      int n_border_vtx = part_ext->vtx_vtx_extended_idx[shift_part+i_part][n_vtx];
      PDM_g_num_t *vtx_ln_to_gn = part_ext->border_vtx_ln_to_gn[shift_part+i_part];
      _shift_ln_to_gn(n_border_vtx, vtx_ln_to_gn, part_ext->shift_by_domain_vtx[i_domain], -1);

      int n_border_face_group          = part_ext->parts[i_domain][i_part].n_face_group;
      PDM_g_num_t *face_group_ln_to_gn = part_ext->border_face_group_ln_to_gn[shift_part+i_part];
      int n_connect_face_group         = part_ext->border_face_group_idx     [shift_part+i_part][n_border_face_group];
      _shift_ln_to_gn(n_connect_face_group, face_group_ln_to_gn, part_ext->shift_by_domain_face_group[i_domain], -1);

    }
    shift_part += part_ext->n_part[i_domain];
  }

}


// static inline
// int
// _is_same_triplet
// (
// int iproc1, int ipart1, int ielt1,
// int iproc2, int ipart2, int ielt2
// )
// {
//   if(iproc1 == iproc2){
//     if(ipart1 == ipart2){
//       if(ielt1 == ielt2){
//         return 1;
//       }
//     }
//   }
//   return 0;
// }

inline
static
int
_lexicographic_equal_int
(
  const int *x,
  const int *y,
  const int          stride
)
{
  int res = x[0] == y[0];
  if(res == 1 && stride > 1) {
    return _lexicographic_equal_int(&x[1], &y[1], stride-1);
  }
  return x[0] == y[0];
}

// static
// int
// _setup_unique_order_triplet
// (
// const int    n_entity,
// const int   *unique_neighbor_entity_idx,
// const int   *unique_neighbor_entity,
//       int  **unique_order_neighbor_entity
// )
// {
//   int* order        = malloc( unique_neighbor_entity_idx[n_entity] * sizeof(int));
//   int* unique_order = malloc( unique_neighbor_entity_idx[n_entity] * sizeof(int));

//   PDM_order_lnum_s(unique_neighbor_entity,
//                    3,
//                    order,
//                    unique_neighbor_entity_idx[n_entity]);

//   int idx_unique = -1;
//   int last_proc  = -1;
//   int last_part  = -1;
//   int last_elmt  = -1;

//   for(int i = 0; i < unique_neighbor_entity_idx[n_entity]; i++){

//     int old_order = order[i];
//     int curr_proc = unique_neighbor_entity[3*old_order  ];
//     int curr_part = unique_neighbor_entity[3*old_order+1];
//     int curr_elmt = unique_neighbor_entity[3*old_order+2];
//     int is_same = _is_same_triplet(last_proc, last_part, last_elmt,
//                                    curr_proc, curr_part, curr_elmt);
//     // log_trace(" curr:: ( %d / %d / %d ) | last:: ( %d / %d / %d ) \n",
//     //             curr_proc, curr_part, curr_elmt,
//     //             last_proc, last_part, last_elmt);

//     if(is_same == 0){ // N'est pas le meme
//       idx_unique++;
//       last_proc = curr_proc;
//       last_part = curr_part;
//       last_elmt = curr_elmt;
//     }
//     unique_order[old_order] = idx_unique;
//   }

//   *unique_order_neighbor_entity = unique_order;
//   free(order);
//   return idx_unique+1;
// }


// static
// void
// _unique_triplet
// (
//   int   n_entity,
//   int  *neighbor_entity_idx,
//   int  *neighbor_entity,
//   int **unique_neighbor_entity_idx,
//   int **unique_neighbor_entity_n,
//   int **unique_neighbor_entity
// )
// {

//   int* _unique_neighbor_entity_idx = malloc( (n_entity + 1) * sizeof(int));
//   int* _unique_neighbor_entity_n   = malloc( (n_entity    ) * sizeof(int));
//   int* _unique_neighbor_entity     = malloc( 3 * neighbor_entity_idx[n_entity] * sizeof(int));
//   int* order                       = malloc(     neighbor_entity_idx[n_entity] * sizeof(int)); // Suralloc

//   _unique_neighbor_entity_idx[0] = 0;
//   for(int i_entity = 0; i_entity < n_entity; ++i_entity) {

//     int beg       = neighbor_entity_idx[i_entity];
//     int n_connect = neighbor_entity_idx[i_entity+1] - beg;

//     PDM_order_lnum_s(&neighbor_entity[3*beg], 3, order, n_connect);

//     _unique_neighbor_entity_n  [i_entity  ] = 0;
//     _unique_neighbor_entity_idx[i_entity+1] = _unique_neighbor_entity_idx[i_entity];

//     int last_proc  = -1;
//     int last_part  = -1;
//     int last_elmt  = -1;
//     for(int i = 0; i < n_connect; ++i) {
//       int old_order   = order[i];
//       int curr_proc   = neighbor_entity[3*(beg+old_order)  ];
//       int curr_part   = neighbor_entity[3*(beg+old_order)+1];
//       int curr_entity = neighbor_entity[3*(beg+old_order)+2];
//       int is_same  = _is_same_triplet(last_proc, last_part, last_elmt,
//                                       curr_proc, curr_part, curr_entity);

//       if(is_same == 0){ // N'est pas le meme
//         // idx_unique++;
//         last_proc = curr_proc;
//         last_part = curr_part;
//         last_elmt = curr_entity;

//         int beg_write = 3 * _unique_neighbor_entity_idx[i_entity+1];
//         // printf("beg_write = %i | curr_proc = %i | curr_part = %i | curr_entity = %i \n", beg_write, curr_proc, curr_part, curr_entity);
//         _unique_neighbor_entity[beg_write  ] = curr_proc;
//         _unique_neighbor_entity[beg_write+1] = curr_part;
//         _unique_neighbor_entity[beg_write+2] = curr_entity;

//         /* Increment the new counter */
//         _unique_neighbor_entity_idx[i_entity+1]++;
//         _unique_neighbor_entity_n  [i_entity  ]++;
//       }
//     }
//   }

//   _unique_neighbor_entity = realloc(_unique_neighbor_entity, 3 * neighbor_entity_idx[n_entity] * sizeof(int));

//   *unique_neighbor_entity_idx = _unique_neighbor_entity_idx;
//   *unique_neighbor_entity_n   = _unique_neighbor_entity_n;
//   *unique_neighbor_entity     = _unique_neighbor_entity;
//   free(order);
// }


static inline
int
_is_same_quadruplet
(
int iproc1, int ipart1, int ielt1, int iinterf1,
int iproc2, int ipart2, int ielt2, int iinterf2
)
{
  if(iproc1 == iproc2){
    if(ipart1 == ipart2){
      if(ielt1 == ielt2){
        if(iinterf1 == iinterf2){
          return 1;
        }
      }
    }
  }
  return 0;
}

static inline
int
_is_same_doublet
(
PDM_g_num_t iproc1, PDM_g_num_t ipart1,
PDM_g_num_t iproc2, PDM_g_num_t ipart2
)
{
  if(iproc1 == iproc2){
    if(ipart1 == ipart2){
      return 1;
    }
  }
  return 0;
}

static
int
_setup_unique_order_quadruplet
(
const int    n_entity,
const int   *unique_neighbor_entity_idx,
const int   *unique_neighbor_entity,
      int  **unique_order_neighbor_entity
)
{
  int* order        = malloc( unique_neighbor_entity_idx[n_entity] * sizeof(int));
  int* unique_order = malloc( unique_neighbor_entity_idx[n_entity] * sizeof(int));

  PDM_order_lnum_s(unique_neighbor_entity,
                   4,
                   order,
                   unique_neighbor_entity_idx[n_entity]);

  int idx_unique = -1;
  int last_proc  = -1;
  int last_part  = -1;
  int last_elmt  = -1;
  int last_inte  = -40;

  for(int i = 0; i < unique_neighbor_entity_idx[n_entity]; i++){

    int old_order = order[i];
    int curr_proc = unique_neighbor_entity[4*old_order  ];
    int curr_part = unique_neighbor_entity[4*old_order+1];
    int curr_elmt = unique_neighbor_entity[4*old_order+2];
    int curr_inte = unique_neighbor_entity[4*old_order+3];
    int is_same = _is_same_quadruplet(last_proc, last_part, last_elmt, last_inte,
                                      curr_proc, curr_part, curr_elmt, curr_inte);
    // log_trace(" curr:: ( %d / %d / %d ) | last:: ( %d / %d / %d ) \n",
    //             curr_proc, curr_part, curr_elmt,
    //             last_proc, last_part, last_elmt);

    if(is_same == 0){ // N'est pas le meme
      idx_unique++;
      last_proc = curr_proc;
      last_part = curr_part;
      last_elmt = curr_elmt;
      last_inte = curr_inte;
    }
    unique_order[old_order] = idx_unique;
  }

  *unique_order_neighbor_entity = unique_order;
  free(order);
  return idx_unique+1;
}


static
void
_unique_quadruplet
(
  int   n_entity,
  int  *neighbor_entity_idx,
  int  *neighbor_entity,
  int **unique_neighbor_entity_idx,
  int **unique_neighbor_entity_n,
  int **unique_neighbor_entity
)
{

  int* _unique_neighbor_entity_idx = malloc( (n_entity + 1) * sizeof(int));
  int* _unique_neighbor_entity_n   = malloc( (n_entity    ) * sizeof(int));
  int* _unique_neighbor_entity     = malloc( 4 * neighbor_entity_idx[n_entity] * sizeof(int));
  int* order                       = malloc(     neighbor_entity_idx[n_entity] * sizeof(int)); // Suralloc

  _unique_neighbor_entity_idx[0] = 0;
  for(int i_entity = 0; i_entity < n_entity; ++i_entity) {

    int beg       = neighbor_entity_idx[i_entity];
    int n_connect = neighbor_entity_idx[i_entity+1] - beg;

    PDM_order_lnum_s(&neighbor_entity[4*beg], 4, order, n_connect);

    _unique_neighbor_entity_n  [i_entity  ] = 0;
    _unique_neighbor_entity_idx[i_entity+1] = _unique_neighbor_entity_idx[i_entity];

    int last_proc  = -1;
    int last_part  = -1;
    int last_elmt  = -1;
    int last_inte  = -40;
    for(int i = 0; i < n_connect; ++i) {
      int old_order   = order[i];
      int curr_proc   = neighbor_entity[4*(beg+old_order)  ];
      int curr_part   = neighbor_entity[4*(beg+old_order)+1];
      int curr_entity = neighbor_entity[4*(beg+old_order)+2];
      int curr_inte   = neighbor_entity[4*(beg+old_order)+3];
      int is_same  = _is_same_quadruplet(last_proc, last_part, last_elmt  , last_inte,
                                         curr_proc, curr_part, curr_entity, curr_inte);

      if(is_same == 0){ // N'est pas le meme
        // idx_unique++;
        last_proc = curr_proc;
        last_part = curr_part;
        last_elmt = curr_entity;
        last_inte = curr_inte;

        int beg_write = 4 * _unique_neighbor_entity_idx[i_entity+1];
        // printf("beg_write = %i | curr_proc = %i | curr_part = %i | curr_entity = %i \n", beg_write, curr_proc, curr_part, curr_entity);
        _unique_neighbor_entity[beg_write  ] = curr_proc;
        _unique_neighbor_entity[beg_write+1] = curr_part;
        _unique_neighbor_entity[beg_write+2] = curr_entity;
        _unique_neighbor_entity[beg_write+3] = curr_inte;

        /* Increment the new counter */
        _unique_neighbor_entity_idx[i_entity+1]++;
        _unique_neighbor_entity_n  [i_entity  ]++;
      }
    }
  }

  _unique_neighbor_entity = realloc(_unique_neighbor_entity, 4 * neighbor_entity_idx[n_entity] * sizeof(int));

  *unique_neighbor_entity_idx = _unique_neighbor_entity_idx;
  *unique_neighbor_entity_n   = _unique_neighbor_entity_n;
  *unique_neighbor_entity     = _unique_neighbor_entity;
  free(order);
}

static
void
triplet_to_quadruplet
(
  int   size,
  int  *triplet,
  int  *array,
  int **quadruplet
)
{
  int *_quadruplet = malloc(4 * size * sizeof(int));
  for(int i = 0; i < size; ++i) {
    _quadruplet[4*i  ] = triplet[3*i  ];
    _quadruplet[4*i+1] = triplet[3*i+1];
    _quadruplet[4*i+2] = triplet[3*i+2];
    _quadruplet[4*i+3] = array  [i];
  }


  *quadruplet = _quadruplet;
}

static
void
quadruplet_to_triplet_and_array
(
  int   size,
  int  *quadruplet,
  int **array,
  int **triplet
)
{
  int *_triplet = malloc(3 * size * sizeof(int));
  int *_array   = malloc(    size * sizeof(int));
  for(int i = 0; i < size; ++i) {
    _triplet[3*i  ] = quadruplet[4*i  ];
    _triplet[3*i+1] = quadruplet[4*i+1];
    _triplet[3*i+2] = quadruplet[4*i+2];
    _array  [i    ] = quadruplet[4*i+3];
  }


  *triplet = _triplet;
  *array   = _array;
}

static
int
_concatenate_neighbor
(
  int             n_part,
  int            *n_entity,
  int           **neighbor_idx,
  int           **init_neighbor_idx,
  int           **init_neighbor_desc,
  int           **next_neighbor_n,
  int           **next_neighbor_desc,
  int          ***concat_neighbor_idx,
  int          ***concat_neighbor_desc
)
{

  int **_concat_neighbor_idx  = malloc(n_part * sizeof(int *));
  int **_concat_neighbor_desc = malloc(n_part * sizeof(int *));

  int is_same = 1;
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int   n_elmt           = n_entity    [i_part];
    int  *_neighbor_idx    = neighbor_idx[i_part];

    int  *_init_neighbor_idx  = init_neighbor_idx [i_part];
    int  *_init_neighbor_desc = init_neighbor_desc[i_part];

    int  *_next_neighbor_n    = next_neighbor_n   [i_part];
    int  *_next_neighbor_desc = next_neighbor_desc[i_part];


    int n_tot_next = 0;
    for(int i = 0; i < n_elmt; ++i) {
      for(int idx = _neighbor_idx[i]; idx < _neighbor_idx[i+1]; ++idx) {
        n_tot_next += _next_neighbor_n[idx];
      }
    }

    int n_concat_idx_tot = _init_neighbor_idx[n_elmt] + n_tot_next;
    _concat_neighbor_idx [i_part] = PDM_array_zeros_int(n_elmt+1);
    _concat_neighbor_desc[i_part] = malloc(4 * n_concat_idx_tot * sizeof(int));

    int idx_read  = 0;
    _concat_neighbor_idx[i_part][0] = 0;
    for(int i = 0; i < n_elmt; ++i) {
      _concat_neighbor_idx[i_part][i+1] = _concat_neighbor_idx[i_part][i];

      /* Current part */
      for(int idx = _init_neighbor_idx[i]; idx < _init_neighbor_idx[i+1]; ++idx) {
        int idx_write = _concat_neighbor_idx[i_part][i+1]++;

        _concat_neighbor_desc[i_part][4*idx_write  ] = _init_neighbor_desc[4*idx  ];
        _concat_neighbor_desc[i_part][4*idx_write+1] = _init_neighbor_desc[4*idx+1];
        _concat_neighbor_desc[i_part][4*idx_write+2] = _init_neighbor_desc[4*idx+2];
        _concat_neighbor_desc[i_part][4*idx_write+3] = _init_neighbor_desc[4*idx+3];

      }

      /* Opp part */
      for(int idx = _neighbor_idx[i]; idx < _neighbor_idx[i+1]; ++idx) {
        for(int k = 0; k < _next_neighbor_n[idx]; ++k) {
          int idx_write = _concat_neighbor_idx[i_part][i+1]++;

          _concat_neighbor_desc[i_part][4*idx_write  ] = _next_neighbor_desc[4*idx_read  ];
          _concat_neighbor_desc[i_part][4*idx_write+1] = _next_neighbor_desc[4*idx_read+1];
          _concat_neighbor_desc[i_part][4*idx_write+2] = _next_neighbor_desc[4*idx_read+2];
          _concat_neighbor_desc[i_part][4*idx_write+3] = _next_neighbor_desc[4*idx_read+3];
          idx_read++;
        }
      }
    }

    /* Unique */
    int *_unique_concat_neighbor_idx  = NULL;
    int *_unique_concat_neighbor_n    = NULL;
    int *_unique_concat_neighbor_desc = NULL;
    _unique_quadruplet(n_elmt,
                       _concat_neighbor_idx [i_part],
                       _concat_neighbor_desc[i_part],
                       &_unique_concat_neighbor_idx,
                       &_unique_concat_neighbor_n,
                       &_unique_concat_neighbor_desc);


    for(int i = 0; i < n_elmt; ++i) {
      int n_old_connect = _init_neighbor_idx[i+1] - _init_neighbor_idx[i];
      if(n_old_connect != _unique_concat_neighbor_n[i]) {
        is_same = 0;
        break;
      }
    }

    // printf("is_same = %i\n", is_same);

    free(_concat_neighbor_idx [i_part]);
    free(_concat_neighbor_desc[i_part]);
    free(_unique_concat_neighbor_n);
    _concat_neighbor_idx [i_part] = _unique_concat_neighbor_idx;
    _concat_neighbor_desc[i_part] = _unique_concat_neighbor_desc;

    /* Debug */
    // PDM_log_trace_graph_nuplet_int(_concat_neighbor_idx[i_part], _concat_neighbor_desc[i_part], 4, n_elmt, "_concat_neighbor_desc :");


  }

  *concat_neighbor_idx  = _concat_neighbor_idx;
  *concat_neighbor_desc = _concat_neighbor_desc;

  return is_same;
}


static
void
_reverse_extended_graph
(
 PDM_MPI_Comm   comm,
 int            n_part_tot,
 int           *part_to_domain,
 int           *shift_part_g_idx,
 int           *n_entity,
 int          **entity_entity_extended_idx,
 int          **entity_entity_extended,
 int          **entity_entity_interface,
 int            n_interface,
 int            n_composed_interface,
 int           *composed_interface_idx,
 int           *composed_interface,
 PDM_g_num_t   *composed_ln_to_gn_sorted,
 int          **revert_entity_extended_idx,
 int          **revert_entity_extended
)
{

  PDM_UNUSED(n_composed_interface);
  PDM_UNUSED(composed_interface_idx);
  PDM_UNUSED(composed_interface);
  PDM_UNUSED(composed_ln_to_gn_sorted);
  PDM_UNUSED(revert_entity_extended_idx);
  PDM_UNUSED(revert_entity_extended);


  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int *send_n = PDM_array_zeros_int(n_rank);

  /* Count */
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    for(int i = 0; i < n_entity[i_part]; ++i) {
      for(int idx = entity_entity_extended_idx[i_part][i]; idx < entity_entity_extended_idx[i_part][i+1]; ++idx) {
        int t_rank = entity_entity_extended[i_part][3*idx];

        int i_interface = entity_entity_interface[i_part][idx]-1;

        if(i_interface < n_interface) {
          send_n[t_rank]++;
        } else { // Composed
          send_n[t_rank]++;
          // int l_interf = PDM_binary_search_long(i_interface+1, composed_ln_to_gn_sorted, n_composed_interface);
          // for(int idx_comp = composed_interface_idx[l_interf]; idx_comp < composed_interface_idx[l_interf+1]; ++idx_comp) {
          //   send_n[t_rank]++;
          // }
        }

      }
    }
  }

  int *send_idx = malloc( (n_rank+1) * sizeof(int));
  send_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx[i+1] = send_idx[i] + send_n[i];
    send_n[i] = 0;
  }

  int *send_buffer = malloc(5 * send_idx[n_rank] * sizeof(int)); // i_part_cur, i_entity_cur, t_part, t_entity, t_interf

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    int i_domain     = part_to_domain[i_part];
    int shift_part_g = shift_part_g_idx[i_domain];
    for(int i = 0; i < n_entity[i_part]; ++i) {
      for(int idx = entity_entity_extended_idx[i_part][i]; idx < entity_entity_extended_idx[i_part][i+1]; ++idx) {
        int t_rank = entity_entity_extended[i_part][3*idx];

        int i_interface = entity_entity_interface[i_part][idx]-1;

        if(i_interface < n_interface) {
          int idx_write = send_idx[t_rank] + send_n[t_rank]++;

          send_buffer[5*idx_write  ] = i_part + shift_part_g;
          send_buffer[5*idx_write+1] = n_entity[i_part] + idx; // Hack here cause we stack at the end the extended ones
          send_buffer[5*idx_write+2] = entity_entity_extended [i_part][3*idx+1];
          send_buffer[5*idx_write+3] = entity_entity_extended [i_part][3*idx+2];
          send_buffer[5*idx_write+4] = entity_entity_interface[i_part][idx];
        } else { // Composed

          int idx_write = send_idx[t_rank] + send_n[t_rank]++;
          send_buffer[5*idx_write  ] = i_part + shift_part_g;
          send_buffer[5*idx_write+1] = n_entity[i_part] + idx; // Hack here cause we stack at the end the extended ones
          send_buffer[5*idx_write+2] = entity_entity_extended [i_part][3*idx+1];
          send_buffer[5*idx_write+3] = entity_entity_extended [i_part][3*idx+2];
          send_buffer[5*idx_write+4] = entity_entity_interface[i_part][idx];
          // int l_interf = PDM_binary_search_long(i_interface+1, composed_ln_to_gn_sorted, n_composed_interface);
          // for(int idx_comp = composed_interface_idx[l_interf]; idx_comp < composed_interface_idx[l_interf+1]; ++idx_comp) {

          //   int idx_write = send_idx[t_rank] + send_n[t_rank]++;
          //   int i_tr = composed_interface[idx_comp];
          //   send_buffer[5*idx_write  ] = i_part + shift_part_g;
          //   send_buffer[5*idx_write+1] = n_entity[i_part] + idx; // Hack here cause we stack at the end the extended ones
          //   send_buffer[5*idx_write+2] = entity_entity_extended[i_part][3*idx+1];
          //   send_buffer[5*idx_write+3] = entity_entity_extended[i_part][3*idx+2];
          //   send_buffer[5*idx_write+4] = i_tr;

          // }
        }

      }
    }
  }

  int *recv_n   = malloc( (n_rank  ) * sizeof(int));
  int *recv_idx = malloc( (n_rank+1) * sizeof(int));
  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT, recv_n, 1, PDM_MPI_INT, comm);


  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }

  int *recv_buffer = malloc(5 * recv_idx[n_rank] * sizeof(int)); // i_part_cur, i_entity_cur, t_part, t_entity

  for(int i = 0; i < n_rank; ++i) {
    send_n  [i] *= 5;
    send_idx[i] *= 5;
    recv_n  [i] *= 5;
    recv_idx[i] *= 5;
  }

  PDM_MPI_Alltoallv (send_buffer,
                     send_n,
                     send_idx,
                     PDM_MPI_INT,
                     recv_buffer,
                     recv_n,
                     recv_idx,
                     PDM_MPI_INT,
                     comm);

  free(send_n);
  free(send_idx);
  free(send_buffer);

  for(int i = 0; i < n_rank; ++i) {
    recv_n  [i] /= 5;
    recv_idx[i] /= 5;
  }
  free(recv_n);

  *revert_entity_extended_idx = recv_idx;
  *revert_entity_extended     = recv_buffer;

}


static
void
_exchange_and_sort_neighbor
(
  PDM_MPI_Comm    comm,
  int             n_part,
  int            *n_entity,
  int           **neighbor_idx,
  int           **neighbor_desc,
  int           **neighbor_interface,
  int           **init_neighbor_idx,
  int           **init_neighbor_desc,
  int          ***all_neighbor_idx,
  int          ***all_neighbor_desc
)
{

  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(comm,
                                                           n_part,
                                                           n_entity,
                                                           neighbor_idx,
                                                           neighbor_desc);

  int **prev_neighbor_idx  = init_neighbor_idx;
  int **prev_neighbor_desc = init_neighbor_desc;

  int is_same    = 0;
  int i_step     = 0;
  int first_step = 1;
  int **concat_neighbor_opp_idx = NULL;
  int **concat_neighbor_opp     = NULL;
  while(is_same != 1) {

    /*
     * Conpute stride
     */
    int **prev_neighbor_n = malloc(n_part * sizeof(int *));
    for(int i_part = 0; i_part < n_part; ++i_part) {
      prev_neighbor_n[i_part] = malloc(n_entity[i_part] * sizeof(int));
      for(int i = 0; i < n_entity[i_part]; ++i) {
        prev_neighbor_n[i_part][i] = prev_neighbor_idx[i_part][i+1] - prev_neighbor_idx[i_part][i];
      }
    }

    int **next_neighbor_opp_n = NULL;
    int **next_neighbor_opp   = NULL;
    PDM_distant_neighbor_exch(dn,
                              4 * sizeof(int),
                              PDM_STRIDE_VAR_INTERLACED,
                              -1,
                              prev_neighbor_n,
                    (void **) prev_neighbor_desc,
                             &next_neighbor_opp_n,
                   (void ***)&next_neighbor_opp);


    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(prev_neighbor_n[i_part]);
    }
    free(prev_neighbor_n);

    is_same = _concatenate_neighbor(n_part,
                                    n_entity,
                                    neighbor_idx,
                                    prev_neighbor_idx,
                                    prev_neighbor_desc,
                                    next_neighbor_opp_n,
                                    next_neighbor_opp,
                                    &concat_neighbor_opp_idx,
                                    &concat_neighbor_opp);

    int is_same_l = is_same;
    PDM_MPI_Allreduce (&is_same_l, &is_same, 1, PDM_MPI_INT, PDM_MPI_MIN, comm);

    /* Init next step */
    if(!first_step) {
      for(int i_part = 0; i_part < n_part; ++i_part) {
        free(prev_neighbor_idx  [i_part]);
        free(prev_neighbor_desc [i_part]);
      }
      free(prev_neighbor_idx);
      free(prev_neighbor_desc);
    }
    prev_neighbor_idx  = concat_neighbor_opp_idx;
    prev_neighbor_desc = concat_neighbor_opp;


    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(next_neighbor_opp_n[i_part]);
      free(next_neighbor_opp  [i_part]);
    }
    free(next_neighbor_opp_n);
    free(next_neighbor_opp);

    first_step = 0;
    i_step++;
    if(i_step > 50) {
      abort();
    }
  }

  PDM_distant_neighbor_free(dn);

  // *all_neighbor_idx  = concat_neighbor_opp_idx;
  // *all_neighbor_desc = concat_neighbor_opp;
  // return;

  /*
   * Filter by removing self and already direct neighbor
   */
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int **filter_neighbor_idx  = malloc(n_part * sizeof(int *));
  int **filter_neighbor_desc = malloc(n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int *_concat_neighbor_opp_idx = concat_neighbor_opp_idx[i_part];
    int *_concat_neighbor_opp     = concat_neighbor_opp    [i_part];

    int *_neighbor_idx   = neighbor_idx      [i_part];
    int *_neighbor_desc  = neighbor_desc     [i_part];
    int *_neighbor_intrf = neighbor_interface[i_part];

    filter_neighbor_idx [i_part] = malloc(     (n_entity[i_part] + 1)                       * sizeof(int));
    filter_neighbor_desc[i_part] = malloc( 4 * (_concat_neighbor_opp_idx[n_entity[i_part]]) * sizeof(int));

    int *_filter_neighbor_idx = filter_neighbor_idx [i_part];
    int *_filter_neighbor     = filter_neighbor_desc[i_part];

    _filter_neighbor_idx[0] = 0;
    for(int i = 0; i < n_entity[i_part]; ++i) {

      _filter_neighbor_idx[i+1] = _filter_neighbor_idx[i];

      for(int idx = _concat_neighbor_opp_idx[i]; idx < _concat_neighbor_opp_idx[i+1]; ++idx) {

        int curr_proc = _concat_neighbor_opp[4*idx  ];
        int curr_part = _concat_neighbor_opp[4*idx+1];
        int curr_elmt = _concat_neighbor_opp[4*idx+2];
        int curr_inte = _concat_neighbor_opp[4*idx+3];

        /* Rm if current elemt is inside */
        int is_define_in_direct_neight = 0;

        if(curr_proc == i_rank && curr_part == i_part && curr_elmt == i) {
          continue;
        }


        /* Brut force */
        for(int idx2 = _neighbor_idx[i]; idx2 < _neighbor_idx[i+1]; ++idx2) {
          int opp_proc = _neighbor_desc[3*idx2  ];
          int opp_part = _neighbor_desc[3*idx2+1];
          int opp_elmt = _neighbor_desc[3*idx2+2];
          int opp_inte = curr_inte; // Normal because i can come from another interface

          is_define_in_direct_neight = _is_same_quadruplet(curr_proc, curr_part, curr_elmt, curr_inte,
                                                           opp_proc , opp_part , opp_elmt , opp_inte);

          if(is_define_in_direct_neight) {
            break;
          }
        }

        if(is_define_in_direct_neight) {
          continue;
        }

        int idx_write = _filter_neighbor_idx[i+1]++;
        _filter_neighbor[4*idx_write  ] = _concat_neighbor_opp[4*idx  ];
        _filter_neighbor[4*idx_write+1] = _concat_neighbor_opp[4*idx+1];
        _filter_neighbor[4*idx_write+2] = _concat_neighbor_opp[4*idx+2];
        _filter_neighbor[4*idx_write+3] = _concat_neighbor_opp[4*idx+3];

      }

      /* On re-rajoute les neighbor */
      for(int idx = _neighbor_idx[i]; idx < _neighbor_idx[i+1]; ++idx) {
        int idx_write = _filter_neighbor_idx[i+1]++;
        _filter_neighbor[4*idx_write  ] = _neighbor_desc[3*idx  ];
        _filter_neighbor[4*idx_write+1] = _neighbor_desc[3*idx+1];
        _filter_neighbor[4*idx_write+2] = _neighbor_desc[3*idx+2];
        _filter_neighbor[4*idx_write+3] = _neighbor_intrf[idx];
      }
    }

    /*
     * Realloc
     */
    filter_neighbor_desc[i_part] = realloc(filter_neighbor_desc[i_part], 4 * (_filter_neighbor_idx[n_entity[i_part]]) * sizeof(int));
    free(_concat_neighbor_opp_idx);
    free(_concat_neighbor_opp);

    if(0 == 1) {
      PDM_log_trace_graph_nuplet_int(_filter_neighbor_idx, filter_neighbor_desc[i_part], 4, n_entity[i_part], "filter_neighbor_desc :");
    }
  }

  free(concat_neighbor_opp_idx);
  free(concat_neighbor_opp    );

  *all_neighbor_idx  = filter_neighbor_idx;
  *all_neighbor_desc = filter_neighbor_desc;
}


static
void
_generate_graph_comm_with_extended
(
 PDM_MPI_Comm                  comm,
 int                           n_part_all_domain,
 int                          *part_to_domain,
 PDM_part_domain_interface_t  *dom_interf,
 PDM_bound_type_t              interface_kind,
 int                          *shift_part_idx,
 int                          *shift_part_g_idx,
 int                          *n_entity,
 int                         **entity_entity_extended_idx,
 int                         **entity_entity_extended,
 int                         **entity_entity_interface,
 int                           n_interface,
 int                           n_composed_interface,
 int                          *composed_interface_idx,
 int                          *composed_interface,
 PDM_g_num_t                  *composed_interface_ln_to_gn,
 PDM_g_num_t                 **border_entity_ln_to_gn,
 int                        ***old_to_new_indices
)
{
  PDM_UNUSED(shift_part_idx);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int *n_entity_extended = malloc(n_part_all_domain * sizeof(int));
  int *n_entity_tot      = malloc(n_part_all_domain * sizeof(int));
  (*old_to_new_indices)  = malloc(n_part_all_domain * sizeof(int *));

  int **neighbor_n         = (int **) malloc( n_part_all_domain * sizeof(int *) );
  int **neighbor_idx       = (int **) malloc( n_part_all_domain * sizeof(int *) );
  int **neighbor_desc      = (int **) malloc( n_part_all_domain * sizeof(int *) );
  int **neighbor_interface = (int **) malloc( n_part_all_domain * sizeof(int *) );

  int **neighbor_opp_n    = (int **) malloc( n_part_all_domain * sizeof(int *) );
  int **neighbor_opp_idx  = (int **) malloc( n_part_all_domain * sizeof(int *) );
  int **neighbor_opp_desc = (int **) malloc( n_part_all_domain * sizeof(int *) );

  for(int i_part = 0; i_part < n_part_all_domain; ++i_part) {
    n_entity_extended[i_part] = entity_entity_extended_idx[i_part][n_entity[i_part]];
    n_entity_tot     [i_part] = n_entity[i_part] + n_entity_extended[i_part];

    neighbor_idx[i_part] = malloc((n_entity_tot[i_part]+1) * sizeof(int));
    neighbor_n  [i_part] = malloc( n_entity_tot[i_part]    * sizeof(int));

    neighbor_opp_idx[i_part] = malloc((n_entity_tot[i_part]+1) * sizeof(int));
    neighbor_opp_n  [i_part] = malloc( n_entity_tot[i_part]    * sizeof(int));

    for(int i = 0; i < n_entity_tot[i_part]; ++i) {
      neighbor_n    [i_part][i] = 0;
      neighbor_opp_n[i_part][i] = 0;
    }
  }

  if(0 == 1) {
    PDM_log_trace_connectivity_int(composed_interface_idx, composed_interface, n_composed_interface, "composed_interface ::");
    PDM_log_trace_array_long(composed_interface_ln_to_gn, n_composed_interface, "composed_interface_ln_to_gn ::");
  }

  int *revert_entity_extended_idx = NULL; // n_rank+1
  int *revert_entity_extended     = NULL; // i_part_cur, i_entity_cur, t_part, t_entity, t_interf
  _reverse_extended_graph(comm,
                          n_part_all_domain,
                          part_to_domain,
                          shift_part_g_idx,
                          n_entity,
                          entity_entity_extended_idx,
                          entity_entity_extended,
                          entity_entity_interface,
                          n_interface,
                          n_composed_interface,
                          composed_interface_idx,
                          composed_interface,
                          composed_interface_ln_to_gn,
                          &revert_entity_extended_idx,
                          &revert_entity_extended);

  if(0 == 1) {
    PDM_log_trace_array_int(revert_entity_extended, 5 * revert_entity_extended_idx[n_rank], "revert_entity_extended :: ");
  }

  /*
   * Create the complet graph with interface and revert extended
   */
  int n_recv_tot = revert_entity_extended_idx[n_rank];
  for(int i = 0; i < n_recv_tot; ++i) {
    int i_part   = revert_entity_extended[5*i+2];
    int i_entity = revert_entity_extended[5*i+3];
    neighbor_n    [i_part][i_entity] += 1;
    neighbor_opp_n[i_part][i_entity] += 1;
  }

  // for(int i_part = 0; i_part < n_part_all_domain; ++i_part) {
  //   for(int i = 0; i < n_entity[i_part]; ++i) {
  //     for(int idx = entity_entity_extended_idx[i_part][i]; idx < entity_entity_extended_idx[i_part][i+1]; ++idx) {
  //       int t_entity = n_entity[i_part] + idx;

  //       int i_interface = entity_entity_interface[i_part][idx]-1;

  //       if(i_interface < n_interface) {
  //         neighbor_n    [i_part][t_entity] += 1;
  //         neighbor_opp_n[i_part][t_entity] += 1;
  //       } else {
  //         neighbor_n    [i_part][t_entity] += 1;
  //         neighbor_opp_n[i_part][t_entity] += 1;
  //         // int l_interf = PDM_binary_search_long(i_interface+1, composed_interface_ln_to_gn, n_composed_interface);
  //         // for(int idx_comp = composed_interface_idx[l_interf]; idx_comp < composed_interface_idx[l_interf+1]; ++idx_comp) {
  //         //   neighbor_n    [i_part][t_entity] += 1;
  //         //   neighbor_opp_n[i_part][t_entity] += 1;
  //         // }
  //       }
  //     }
  //   }
  // }


  /*
   * Loop over all interfaces to create distant neighbor structure
   */
  for(int i_part = 0; i_part < n_part_all_domain; ++i_part) {

    int i_domain     = part_to_domain[i_part];
    int shift_part_g = shift_part_g_idx[i_domain];

    int* _neighbor_n   = neighbor_n    [i_part];
    int* _neighbor_idx = neighbor_idx  [i_part];

    int* _neighbor_opp_n   = neighbor_opp_n    [i_part];
    int* _neighbor_opp_idx = neighbor_opp_idx  [i_part];

    int           *interface_pn       = malloc(n_interface * sizeof(int          ));
    PDM_g_num_t  **interface_ln_to_gn = malloc(n_interface * sizeof(PDM_g_num_t *));
    int          **interface_sgn      = malloc(n_interface * sizeof(int         *));
    int          **interface_sens     = malloc(n_interface * sizeof(int         *));
    int          **interface_ids      = malloc(n_interface * sizeof(int         *));
    int          **interface_ids_idx  = malloc(n_interface * sizeof(int         *));
    int          **interface_dom      = malloc(n_interface * sizeof(int         *));
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      PDM_part_domain_interface_get(dom_interf,
                                    interface_kind,
                                    i_domain,
                                    i_part,
                                    i_interface,
                                    &interface_pn      [i_interface],
                                    &interface_ln_to_gn[i_interface],
                                    &interface_sgn     [i_interface],
                                    &interface_sens    [i_interface],
                                    &interface_ids     [i_interface],
                                    &interface_ids_idx [i_interface],
                                    &interface_dom     [i_interface]);
    }

    /*
     * First step : Count interface to add in distant neighbor due to connectivity betwenn domain
     */
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

      // log_trace("-------------------------------- i_interface = %i  -------------------------------- \n", i_interface);
      // PDM_log_trace_array_int(interface_sgn[i_interface], interface_pn[i_interface], "interface_sgn :: ");

      for(int idx_entity = 0; idx_entity < interface_pn[i_interface]; ++idx_entity) {

        // Search the first in list that is in current part/proc
        // int i_proc_cur   = -1;
        // int i_part_cur   = -1;
        int i_entity_cur = -1;
        int found        = 0;
        int idx_current  = -1;
        for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
          int i_proc_opp   = interface_ids[i_interface][3*j  ];
          int i_part_opp   = interface_ids[i_interface][3*j+1];
          int i_entity_opp = interface_ids[i_interface][3*j+2];

          if(i_proc_opp == i_rank && i_part_opp == i_part && found == 0) {
            // i_proc_cur   = i_proc_opp;
            // i_part_cur   = i_part_opp;
            i_entity_cur = i_entity_opp;
            idx_current  = j;
            found = 1;
          }
        }

        if(!found) {
          continue;
        }

        // Il manque une notion de direction sinon on sait pas dans quelle sens va le raccord

        assert(found == 1);

        // log_trace("i_proc_cur = %i | i_part_cur = %i | i_entity_cur = %i | sgn = %i \n", i_proc_cur, i_part_cur, i_entity_cur, interface_sgn[i_interface][idx_entity]);

        // Only add the opposite part of the graph
        for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
          // int i_proc_opp   = interface_ids[i_interface][3*j  ];
          // int i_part_opp   = interface_ids[i_interface][3*j+1];
          // int i_entity_opp = interface_ids[i_interface][3*j+2];

          if(idx_current != j) {
            // log_trace("\t i_proc_opp = %i | i_part_opp = %i | i_entity_opp = %i \n", i_proc_opp, i_part_opp, i_entity_opp);
            _neighbor_n[i_entity_cur] += 1;
            _neighbor_opp_n[i_entity_cur] += 1;
          // } else {
          //   _neighbor_opp_n[i_entity_cur] += 1;
          }
        }
      }
    }

    /*
     * Rajout des connexions entre partitions !!!
     */

    /* Compute index */
    _neighbor_idx[0] = 0;
    _neighbor_opp_idx[0] = 0;
    for(int i_entity = 0; i_entity < n_entity_tot[i_part]; ++i_entity) {
      _neighbor_idx    [i_entity+1] = _neighbor_idx    [i_entity] + _neighbor_n    [i_entity];
      _neighbor_opp_idx[i_entity+1] = _neighbor_opp_idx[i_entity] + _neighbor_opp_n[i_entity];
      _neighbor_n    [i_entity] = 0;
      _neighbor_opp_n[i_entity] = 0;
    }

    neighbor_desc     [i_part] = (int *) malloc( 3 * _neighbor_idx    [n_entity_tot[i_part]] * sizeof(int) );
    neighbor_opp_desc [i_part] = (int *) malloc( 4 * _neighbor_opp_idx[n_entity_tot[i_part]] * sizeof(int) );
    neighbor_interface[i_part] = (int *) malloc(     _neighbor_idx[n_entity_tot[i_part]] * sizeof(int) );
    int* _neighbor_desc      = neighbor_desc     [i_part];
    int* _neighbor_opp_desc  = neighbor_opp_desc [i_part];
    int* _neighbor_interface = neighbor_interface[i_part];

    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      for(int idx_entity = 0; idx_entity < interface_pn[i_interface]; ++idx_entity) {

        // Search the first in list that is in current part/proc
        // int i_proc_cur   = -1;
        // int i_part_cur   = -1;
        int i_entity_cur = -1;
        int found        = 0;
        int idx_current  = -1;
        for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
          int i_proc_opp   = interface_ids[i_interface][3*j  ];
          int i_part_opp   = interface_ids[i_interface][3*j+1];
          int i_entity_opp = interface_ids[i_interface][3*j+2];

          if(i_proc_opp == i_rank && i_part_opp == i_part && found == 0) {
            // i_proc_cur   = i_proc_opp;
            // i_part_cur   = i_part_opp;
            i_entity_cur = i_entity_opp;
            idx_current  = j;
            found = 1;
          }
        }

        if(!found) {
          continue;
        }

        assert(found == 1);

        // Only add the opposite part of the graph
        for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
          int i_proc_opp   = interface_ids[i_interface][3*j  ];
          int i_part_opp   = interface_ids[i_interface][3*j+1];
          int i_entity_opp = interface_ids[i_interface][3*j+2];

          if(idx_current != j) {
            // log_trace("\t i_proc_opp = %i | i_part_opp = %i | i_entity_opp = %i \n", i_proc_opp, i_part_opp, i_entity_opp);
            int idx_write = _neighbor_idx[i_entity_cur] + _neighbor_n[i_entity_cur]++;
            _neighbor_desc[3*idx_write  ] = i_proc_opp;              // i_proc_opp;
            _neighbor_desc[3*idx_write+1] = i_part_opp+shift_part_g; // i_part_opp
            _neighbor_desc[3*idx_write+2] = i_entity_opp;            // i_entity_opp
            _neighbor_interface[idx_write] = (i_interface+1) * interface_sgn[i_interface][idx_entity];

            idx_write = _neighbor_opp_idx[i_entity_cur] + _neighbor_opp_n[i_entity_cur]++;
            _neighbor_opp_desc[4*idx_write  ] = i_proc_opp;              // i_proc_opp;
            _neighbor_opp_desc[4*idx_write+1] = i_part_opp+shift_part_g; // i_part_opp
            _neighbor_opp_desc[4*idx_write+2] = i_entity_opp;            // i_entity_opp
            _neighbor_opp_desc[4*idx_write+3] = (i_interface+1) * interface_sgn[i_interface][idx_entity];
          }
        }
      }

      free(interface_pn      );
      free(interface_ln_to_gn);
      free(interface_sgn     );
      free(interface_sens    );
      free(interface_ids     );
      free(interface_ids_idx );
      free(interface_dom     );
    }
  }

  /* Remplissage de fin */

  // for(int i_part = 0; i_part < n_part_all_domain; ++i_part) {
  //   for(int i = 0; i < n_entity[i_part]; ++i) {
  //     for(int idx = entity_entity_extended_idx[i_part][i]; idx < entity_entity_extended_idx[i_part][i+1]; ++idx) {
  //       int t_entity = n_entity[i_part] + idx;

  //       int i_interface = entity_entity_interface[i_part][idx]-1;

  //       if(i_interface < n_interface) {

  //         int idx_write     = neighbor_idx    [i_part][t_entity] + neighbor_n    [i_part][t_entity]++;
  //         int idx_write_opp = neighbor_opp_idx[i_part][t_entity] + neighbor_opp_n[i_part][t_entity]++;

  //         neighbor_desc     [i_part][3*idx_write  ] = entity_entity_extended[i_part][3*idx  ];
  //         neighbor_desc     [i_part][3*idx_write+1] = entity_entity_extended[i_part][3*idx+1];
  //         neighbor_desc     [i_part][3*idx_write+2] = entity_entity_extended[i_part][3*idx+2];
  //         neighbor_interface[i_part][idx_write    ] = entity_entity_interface[i_part][idx];

  //         neighbor_opp_desc[i_part][4*idx_write_opp  ] = entity_entity_extended [i_part][3*idx  ];
  //         neighbor_opp_desc[i_part][4*idx_write_opp+1] = entity_entity_extended [i_part][3*idx+1];
  //         neighbor_opp_desc[i_part][4*idx_write_opp+2] = entity_entity_extended [i_part][3*idx+2];
  //         neighbor_opp_desc[i_part][4*idx_write_opp+3] = entity_entity_interface[i_part][idx];
  //       } else {

  //         int idx_write     = neighbor_idx    [i_part][t_entity] + neighbor_n    [i_part][t_entity]++;
  //         int idx_write_opp = neighbor_opp_idx[i_part][t_entity] + neighbor_opp_n[i_part][t_entity]++;
  //         neighbor_desc     [i_part][3*idx_write  ] = entity_entity_extended[i_part][3*idx  ];
  //         neighbor_desc     [i_part][3*idx_write+1] = entity_entity_extended[i_part][3*idx+1];
  //         neighbor_desc     [i_part][3*idx_write+2] = entity_entity_extended[i_part][3*idx+2];
  //         neighbor_interface[i_part][idx_write    ] = entity_entity_interface[i_part][idx];

  //         neighbor_opp_desc[i_part][4*idx_write_opp  ] = entity_entity_extended [i_part][3*idx  ];
  //         neighbor_opp_desc[i_part][4*idx_write_opp+1] = entity_entity_extended [i_part][3*idx+1];
  //         neighbor_opp_desc[i_part][4*idx_write_opp+2] = entity_entity_extended [i_part][3*idx+2];
  //         neighbor_opp_desc[i_part][4*idx_write_opp+3] = entity_entity_interface[i_part][idx];
  //         // int l_interf = PDM_binary_search_long(i_interface+1, composed_interface_ln_to_gn, n_composed_interface);
  //         // for(int idx_comp = composed_interface_idx[l_interf]; idx_comp < composed_interface_idx[l_interf+1]; ++idx_comp) {

  //         //   int idx_write     = neighbor_idx    [i_part][t_entity] + neighbor_n    [i_part][t_entity]++;
  //         //   int idx_write_opp = neighbor_opp_idx[i_part][t_entity] + neighbor_opp_n[i_part][t_entity]++;
  //         //   int i_tr = composed_interface[idx_comp];

  //         //   neighbor_desc     [i_part][3*idx_write  ] = entity_entity_extended[i_part][3*idx  ];
  //         //   neighbor_desc     [i_part][3*idx_write+1] = entity_entity_extended[i_part][3*idx+1];
  //         //   neighbor_desc     [i_part][3*idx_write+2] = entity_entity_extended[i_part][3*idx+2];
  //         //   neighbor_interface[i_part][idx_write    ] = i_tr;

  //         //   neighbor_opp_desc[i_part][4*idx_write_opp  ] = entity_entity_extended [i_part][3*idx  ];
  //         //   neighbor_opp_desc[i_part][4*idx_write_opp+1] = entity_entity_extended [i_part][3*idx+1];
  //         //   neighbor_opp_desc[i_part][4*idx_write_opp+2] = entity_entity_extended [i_part][3*idx+2];
  //         //   neighbor_opp_desc[i_part][4*idx_write_opp+3] = i_tr;
  //         // }
  //       }
  //     }
  //   }
  // }

  for(int t_rank = 0; t_rank < n_rank; ++t_rank) {
    for(int i = revert_entity_extended_idx[t_rank]; i < revert_entity_extended_idx[t_rank+1]; ++i) {
      int i_part   = revert_entity_extended[5*i+2];
      int i_entity = revert_entity_extended[5*i+3];
      int idx_write     = neighbor_idx    [i_part][i_entity] + neighbor_n    [i_part][i_entity]++;
      int idx_write_opp = neighbor_opp_idx[i_part][i_entity] + neighbor_opp_n[i_part][i_entity]++;

      neighbor_desc[i_part][3*idx_write  ] = t_rank;
      neighbor_desc[i_part][3*idx_write+1] = revert_entity_extended[5*i  ];
      neighbor_desc[i_part][3*idx_write+2] = revert_entity_extended[5*i+1];
      neighbor_interface[i_part][idx_write] = -revert_entity_extended[5*i+4];

      neighbor_opp_desc[i_part][4*idx_write_opp  ] = t_rank;
      neighbor_opp_desc[i_part][4*idx_write_opp+1] = revert_entity_extended[5*i  ];
      neighbor_opp_desc[i_part][4*idx_write_opp+2] = revert_entity_extended[5*i+1];
      neighbor_opp_desc[i_part][4*idx_write_opp+3] = -revert_entity_extended[5*i+4];

    }
  }

  free(revert_entity_extended_idx);
  free(revert_entity_extended);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part_all_domain; ++i_part) {
      PDM_log_trace_graph_nuplet_int(neighbor_opp_idx[i_part],  neighbor_opp_desc[i_part], 4, n_entity_tot[i_part], "neighbor_opp_desc :");
    }

  }

  int **neighbor_entity_idx  = NULL;
  int **neighbor_entity_desc = NULL;
  _exchange_and_sort_neighbor(comm,
                              n_part_all_domain,
                              n_entity_tot,
                              neighbor_idx,
                              neighbor_desc,
                              neighbor_interface,
                              neighbor_opp_idx,
                              neighbor_opp_desc,
                              &neighbor_entity_idx,
                              &neighbor_entity_desc);

  for(int i_part = 0; i_part < n_part_all_domain; ++i_part) {
    free(neighbor_opp_n    [i_part]);
    free(neighbor_idx      [i_part]);
    free(neighbor_n        [i_part]);
    free(neighbor_desc     [i_part]);
    free(neighbor_opp_idx  [i_part]);
    free(neighbor_opp_desc [i_part]);
    free(neighbor_interface[i_part]);
  }
  free(neighbor_n);
  free(neighbor_idx);
  free(neighbor_opp_n);
  free(neighbor_desc);
  free(neighbor_opp_idx);
  free(neighbor_opp_desc);
  free(neighbor_interface);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part_all_domain; ++i_part) {
      PDM_log_trace_graph_nuplet_int(neighbor_entity_idx[i_part],  neighbor_entity_desc[i_part], 4, n_entity_tot[i_part], "neighbor_entity_desc :");
    }
  }

  /*
   * Tentative -> On cherche tout le nuplet identique dans les extended
   */
  for(int i_part = 0; i_part < n_part_all_domain; ++i_part) {

    int *sort_neighbor_n = malloc(n_entity_extended[i_part] * sizeof(int));
    int *order           = malloc(n_entity_extended[i_part] * sizeof(int));

    int *_neighbor_entity_idx  = neighbor_entity_idx [i_part];
    int *_neighbor_entity_desc = neighbor_entity_desc[i_part];

    int size_neight = 0;
    for(int i = n_entity[i_part]; i < n_entity_tot[i_part]; ++i) {
      sort_neighbor_n[i-n_entity[i_part]] = _neighbor_entity_idx[i+1] - _neighbor_entity_idx[i];
      order     [i-n_entity[i_part]]      = i; //-n_entity[i_part];
      size_neight                        += _neighbor_entity_idx[i+1] - _neighbor_entity_idx[i];
    }

    PDM_sort_int(sort_neighbor_n, order, n_entity_extended[i_part]);

    int *sort_neighbor_idx = malloc( (n_entity_extended[i_part]+1) * sizeof(int));
    sort_neighbor_idx[0] = 0;
    for(int i = 0; i < n_entity_extended[i_part]; ++i) {
      sort_neighbor_idx[i+1] = sort_neighbor_idx[i] + sort_neighbor_n[i];
    }

    int *sort1_neighbor = malloc(4 * size_neight * sizeof(int));

    for(int i = 0; i < n_entity_extended[i_part]; ++i) {
      int old_order = order[i];

      assert(_neighbor_entity_idx[old_order+1]-_neighbor_entity_idx[old_order] == sort_neighbor_n[i]);
      for(int k = 0; k < sort_neighbor_n[i]; ++k) {
        int idx_read  = _neighbor_entity_idx[old_order] + k ;
        int idx_write = sort_neighbor_idx[i] + k;


        sort1_neighbor[4*idx_write  ] = _neighbor_entity_desc[4*idx_read  ];
        sort1_neighbor[4*idx_write+1] = _neighbor_entity_desc[4*idx_read+1];
        sort1_neighbor[4*idx_write+2] = _neighbor_entity_desc[4*idx_read+2];
        sort1_neighbor[4*idx_write+3] = _neighbor_entity_desc[4*idx_read+3];
      }
    }

    if(0 == 1) {
      PDM_log_trace_graph_nuplet_int(sort_neighbor_idx,  sort1_neighbor, 4, n_entity_extended[i_part], "sort1_neighbor :");
    }

    /* Identify bucket size */
    int n_bucket = -1;
    int *bucket_idx = malloc( (n_entity_extended[i_part]+1) * sizeof(int));
    int first_size = -1;
    for(int i =  0; i < n_entity_extended[i_part]; ++i) {
      if(sort_neighbor_n[i] > first_size) {
        bucket_idx[n_bucket+1] = i;
        first_size = sort_neighbor_n[i];
        n_bucket++;
      }
    }
    bucket_idx[ n_bucket+1] = n_entity_extended[i_part];
    n_bucket++;
    bucket_idx = realloc(bucket_idx,  (n_bucket+1) * sizeof(int));

    // PDM_log_trace_array_int(bucket_idx,  n_bucket+1, "bucket_idx :");
    // PDM_log_trace_array_int(sort_neighbor_n,  n_entity_extended[i_part], "sort_neighbor_n :");

    int max_cst_size = -1;
    for(int i_bucket = 0; i_bucket < n_bucket; ++i_bucket) {
      int beg_bucket = bucket_idx[i_bucket];
      int cst_size   = 4 * sort_neighbor_n[beg_bucket];

      max_cst_size = PDM_MAX(max_cst_size, cst_size);

      int beg = 4 * sort_neighbor_idx[beg_bucket];
      int n_in_bucket = bucket_idx[i_bucket+1] - beg_bucket;

      int *bucket_order = malloc(n_in_bucket * sizeof(int));

      PDM_order_lnum_s(&sort1_neighbor[beg], cst_size, bucket_order, n_in_bucket);
      PDM_order_array (n_in_bucket, cst_size * sizeof(int), bucket_order, &sort1_neighbor[beg]);
      PDM_order_array (n_in_bucket, sizeof(int), bucket_order, &order[beg_bucket]);

      // PDM_log_trace_array_int(bucket_order,  n_in_bucket, "bucket_order :");

      free(bucket_order);

    }


    /*
     * Shift order
     */
    for(int i = 0; i < n_entity_extended[i_part]; ++i) {
      order[i] = order[i] - n_entity[i_part];
    }

    if(0 == 1) {
      PDM_log_trace_graph_nuplet_int(sort_neighbor_idx,  sort1_neighbor, 4, n_entity_extended[i_part], "sort1_neighbor (final) :");

      PDM_log_trace_array_int(order,  n_entity_extended[i_part], "order :");
    }

    // Calcul des idx qui ont le meme type (donc dans un bucket les unifié)
    // Pour chaque sous bucket on unifie si le gnum est différent et on les associes

    int* unique_idx = malloc((n_entity_extended[i_part]+1) * sizeof(int));

    for(int i = 0; i < n_entity_extended[i_part]+1; ++i) {
      unique_idx[i] = 0;
    }

    int* last_value = malloc((max_cst_size               ) * sizeof(int));
    int i_unique = 0;
    for(int i_bucket = 0; i_bucket < n_bucket; ++i_bucket) {

      int beg_bucket = bucket_idx[i_bucket];
      int cst_size   = 4 * sort_neighbor_n[beg_bucket];
      int n_in_bucket = bucket_idx[i_bucket+1] - beg_bucket;

      int beg = 4 * sort_neighbor_idx[beg_bucket];

      for(int j = 0; j < (int) cst_size; ++j) {
        last_value[j] = sort1_neighbor[beg+j];
      }

      for(int k = 0; k < n_in_bucket; ++k) {
        int is_same = _lexicographic_equal_int(last_value, &sort1_neighbor[beg+cst_size*k], cst_size);
        if(is_same == 0){ // N'est pas le meme
          for(int j = 0; j < (int) cst_size; ++j) {
            last_value[j] = sort1_neighbor[beg+cst_size*k+j];
          }

          i_unique++;
          unique_idx[i_unique+1]++;

        } else {
          unique_idx[i_unique+1]++;
        }
      }

      i_unique++;
    }

    for(int i = 0; i < i_unique; ++i) {
      unique_idx[i+1] += unique_idx[i];
    }

    if(1 == 0) {
      PDM_log_trace_array_int(unique_idx,  i_unique+1, "unique_idx :");
      for(int i = 0; i < n_entity_extended[i_part]; ++i) {
        int old_order = order[i];
        printf("[%i] gnum = "PDM_FMT_G_NUM"\n", i, border_entity_ln_to_gn[i_part][old_order]);
      }
    }


    int* is_treated = PDM_array_zeros_int(n_entity_extended[i_part]);

    int* _old_to_new_indices = malloc(n_entity_extended[i_part] * sizeof(int));
    for(int i = 0; i < n_entity_extended[i_part]; ++i) {
      _old_to_new_indices[i] = -10000000;
    }

    int new_indices = 0;
    for(int i = 0; i < i_unique; ++i) {

      // n*n algo but n is supposed to be max at 3
      for(int idx = unique_idx[i]; idx < unique_idx[i+1]; ++idx) {

        if(is_treated[idx] == 1) {
          continue;
        }
        int old_order = order[idx];
        PDM_g_num_t gnum1 = border_entity_ln_to_gn[i_part][old_order];

        for(int idx2 = unique_idx[i]+1; idx2 < unique_idx[i+1]; ++idx2) {
          if(is_treated[idx2] == 1) {
            continue;
          }
          int old_order2 = order[idx2];
          PDM_g_num_t gnum2 = border_entity_ln_to_gn[i_part][old_order2];

          if(gnum1 != gnum2) {
            // printf("Is same !!! --> %i %i | %i %i --> new_indices = %i \n", old_order, old_order2, (int)gnum1, (int)gnum2, new_indices);
            is_treated[idx ] = 1;
            is_treated[idx2] = 1;
            _old_to_new_indices[old_order ] = new_indices;
            _old_to_new_indices[old_order2] = new_indices;
            // break; // This one is false
          }
        }

        if(is_treated[idx ] == 1) {
          new_indices++;
        } else{
          _old_to_new_indices[old_order] = new_indices++;
        }
      }


    }

    // PDM_log_trace_array_int(_old_to_new_indices, n_entity_extended[i_part], "_old_to_new_indices :");

    (*old_to_new_indices[i_part]) = _old_to_new_indices;




    free(is_treated);
    free(unique_idx);
    free(last_value);
    free(bucket_idx);
    free(order);
    free(sort1_neighbor);
    free(sort_neighbor_n);
    free(sort_neighbor_idx);
  }




  for(int i_part = 0; i_part < n_part_all_domain; ++i_part) {
    free(neighbor_entity_idx[i_part]);
    free(neighbor_entity_desc[i_part]);
  }
  free(neighbor_entity_idx);
  free(neighbor_entity_desc);


  free(n_entity_extended);
  free(n_entity_tot);

}


static
void
_create_cell_cell_graph
(
  PDM_part_extension_t *part_ext,
  PDM_extend_type_t     extend_type
)
{
  // assert(extend_type == PDM_EXTEND_FROM_FACE);

  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  part_ext->cell_cell_idx = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->cell_cell     = (int **) malloc( n_part_loc_all_domain * sizeof(int *));

  if(extend_type == PDM_EXTEND_FROM_FACE ){

    /* In order to call generic fonction we tranform face_cell with idx */
    int ***face_cell_idx = malloc( part_ext->n_domain * sizeof(int ***));
    int ***face_cell     = malloc( part_ext->n_domain * sizeof(int ***));

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

      int  *pn_face        = (int * ) malloc( part_ext->n_part[i_domain] * sizeof(int  ));
      int  *pn_cell        = (int * ) malloc( part_ext->n_part[i_domain] * sizeof(int  ));
      int **pface_cell     = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
      int **pcell_face     = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));
      int **pcell_face_idx = (int **) malloc( part_ext->n_part[i_domain] * sizeof(int *));

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_face       [i_part] = part_ext->parts[i_domain][i_part].n_face;
        pn_cell       [i_part] = part_ext->parts[i_domain][i_part].n_cell;
        pface_cell    [i_part] = part_ext->parts[i_domain][i_part].face_cell;
        pcell_face_idx[i_part] = part_ext->parts[i_domain][i_part].cell_face_idx;
        pcell_face    [i_part] = part_ext->parts[i_domain][i_part].cell_face;
      }

      PDM_part_connectivity_to_connectity_idx(part_ext->n_part[i_domain],
                                              pn_face,
                                              pface_cell,
                                              &face_cell_idx[i_domain],
                                              &face_cell[i_domain]);


      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

        // PDM_log_trace_array_int(face_cell_idx[i_domain][i_part], pn_face[i_part], "face_cell_idx::");
        // PDM_log_trace_array_int(face_cell    [i_domain][i_part], face_cell_idx[i_domain][i_part][pn_face[i_part]], "face_cell::");
        // printf("[%i] face_cell_idx -> \n ", i_part);
        // for(int i_face = 0; i_face < pn_face[i_part]; ++i_face) {
        //   printf(" (%i) = ", i_face);
        //   for(int idx_cell = face_cell_idx[i_domain][i_part][i_face]; idx_cell < face_cell_idx[i_domain][i_part][i_face+1]; ++idx_cell) {
        //     printf(" %i", face_cell[i_domain][i_part][idx_cell]-1);
        //   }
        //   printf("\n");
        // }

        // On fait égalemnt le cell_cell
        // printf("PDM_combine_connectivity[%i] -> %i %i \n", i_part, pn_cell[i_part], pn_face[i_part]);
        PDM_combine_connectivity(pn_cell[i_part],
                                 pcell_face_idx[i_part],
                                 pcell_face[i_part],
                                 face_cell_idx[i_domain][i_part],
                                 face_cell[i_domain][i_part],
                                 &part_ext->cell_cell_idx[i_part+shift_part],
                                 &part_ext->cell_cell[i_part+shift_part]);

        // Remove sign for cell_cell
        for(int i = 0; i < part_ext->cell_cell_idx[i_part+shift_part][pn_cell[i_part]]; ++i) {
          part_ext->cell_cell[i_part+shift_part][i] = PDM_ABS(part_ext->cell_cell[i_part+shift_part][i]);
        }

        /*
         * Setup shortcut and free useless memory
         */
        part_ext->entity_cell_idx[i_part+shift_part] = face_cell_idx[i_domain][i_part];
        part_ext->entity_cell    [i_part+shift_part] = face_cell    [i_domain][i_part];

        part_ext->entity_cell_n  [i_part+shift_part] = (int * ) malloc( pn_face[i_part] * sizeof(int));
        for(int i_face = 0; i_face < pn_face[i_part]; ++i_face) {
          int n_connected = part_ext->entity_cell_idx[i_part+shift_part][i_face+1] - part_ext->entity_cell_idx[i_part+shift_part][i_face];
          part_ext->entity_cell_n[i_part+shift_part][i_face] = n_connected;
        }
      }

      free(pn_face);
      free(pn_cell);
      free(pface_cell);
      free(pcell_face);
      free(pcell_face_idx);

      free(face_cell_idx[i_domain]);
      free(face_cell[i_domain]);

      shift_part += part_ext->n_part[i_domain];
    }

    free(face_cell_idx);
    free(face_cell);
  } else if(extend_type == PDM_EXTEND_FROM_VTX) {

    // Check

    // Compute cell_cell = cell_face + face_edge + edge_vtx then transpose

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        int pn_cell = part_ext->parts[i_domain][i_part].n_cell;
        int pn_vtx  = part_ext->parts[i_domain][i_part].n_vtx;

        int* cell_vtx_idx = NULL;
        int* cell_vtx     = NULL;

        if (part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE]) {
          assert(part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
          assert(part_ext->parts[i_domain][i_part].face_edge_idx != NULL);
          assert(part_ext->parts[i_domain][i_part].face_edge     != NULL);

          int pn_edge = part_ext->parts[i_domain][i_part].n_edge;
          int *edge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, pn_edge);

          /* Compute cell_edge */
          int* cell_edge_idx = NULL;
          int* cell_edge     = NULL;
          PDM_combine_connectivity(part_ext->parts[i_domain][i_part].n_cell,
                                   part_ext->parts[i_domain][i_part].cell_face_idx,
                                   part_ext->parts[i_domain][i_part].cell_face,
                                   part_ext->parts[i_domain][i_part].face_edge_idx,
                                   part_ext->parts[i_domain][i_part].face_edge,
                                   &cell_edge_idx,
                                   &cell_edge);

          /* Compute cell_vtx */
          PDM_combine_connectivity(part_ext->parts[i_domain][i_part].n_cell,
                                   cell_edge_idx,
                                   cell_edge,
                                   edge_vtx_idx,
                                   part_ext->parts[i_domain][i_part].edge_vtx,
                                   &cell_vtx_idx,
                                   &cell_vtx);
          free(cell_edge_idx);
          free(cell_edge);
          free(edge_vtx_idx);
        }
        else {
          assert(part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX]);
          assert(part_ext->parts[i_domain][i_part].face_vtx_idx != NULL);

          /* Compute cell_vtx */
          PDM_combine_connectivity(part_ext->parts[i_domain][i_part].n_cell,
                                   part_ext->parts[i_domain][i_part].cell_face_idx,
                                   part_ext->parts[i_domain][i_part].cell_face,
                                   part_ext->parts[i_domain][i_part].face_vtx_idx,
                                   part_ext->parts[i_domain][i_part].face_vtx,
                                   &cell_vtx_idx,
                                   &cell_vtx);
        }


        int *vtx_cell     = NULL;
        int *vtx_cell_idx = NULL;
        PDM_connectivity_transpose(part_ext->parts[i_domain][i_part].n_cell,
                                   part_ext->parts[i_domain][i_part].n_vtx,
                                   cell_vtx_idx,
                                   cell_vtx,
                                   &vtx_cell_idx,
                                   &vtx_cell);

        PDM_combine_connectivity(part_ext->parts[i_domain][i_part].n_cell,
                                 cell_vtx_idx,
                                 cell_vtx,
                                 vtx_cell_idx,
                                 vtx_cell,
                                 &part_ext->cell_cell_idx[i_part+shift_part],
                                 &part_ext->cell_cell[i_part+shift_part]);

        // Remove sign for cell_cell
        for(int i = 0; i < part_ext->cell_cell_idx[i_part+shift_part][pn_cell]; ++i) {
          part_ext->cell_cell[i_part+shift_part][i] = PDM_ABS(part_ext->cell_cell[i_part+shift_part][i]);
        }
        free(cell_vtx_idx);
        free(cell_vtx);

        /*
         * Setup shortcut and free useless memory
         */
        part_ext->entity_cell_idx[i_part+shift_part] = vtx_cell_idx;
        part_ext->entity_cell    [i_part+shift_part] = vtx_cell    ;
        part_ext->entity_cell_n  [i_part+shift_part] = (int * ) malloc( pn_vtx * sizeof(int));
        for(int i_vtx = 0; i_vtx < pn_vtx; ++i_vtx) {
          int n_connected = part_ext->entity_cell_idx[i_part+shift_part][i_vtx+1] - part_ext->entity_cell_idx[i_part+shift_part][i_vtx];
          part_ext->entity_cell_n[i_part+shift_part][i_vtx] = n_connected;
        }

      } /* End i_part */
      shift_part += part_ext->n_part[i_domain];
    }


  } else {
    abort();
  }
}

static
void
_compute_other_domain_interface
(
  PDM_part_extension_t *part_ext
)
{
  if(part_ext->pdi == NULL) {
    return;
  }

  int is_describe_vtx  = PDM_part_domain_interface_exist_get(part_ext->pdi, PDM_BOUND_TYPE_VTX );
  int is_describe_edge = PDM_part_domain_interface_exist_get(part_ext->pdi, PDM_BOUND_TYPE_EDGE);
  int is_describe_face = PDM_part_domain_interface_exist_get(part_ext->pdi, PDM_BOUND_TYPE_FACE);

  int is_describe_vtx_l  = is_describe_vtx;
  int is_describe_edge_l = is_describe_edge;
  int is_describe_face_l = is_describe_face;
  PDM_MPI_Allreduce(&is_describe_vtx_l , &is_describe_vtx , 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);
  PDM_MPI_Allreduce(&is_describe_edge_l, &is_describe_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);
  PDM_MPI_Allreduce(&is_describe_face_l, &is_describe_face, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

  int have_edge = 0;
  int have_face = 0;

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      if(part_ext->parts[i_domain][i_part].n_edge > 0 &&
         part_ext->parts[i_domain][i_part].edge_ln_to_gn != NULL) {
        have_edge = 1;
      }
      if(part_ext->parts[i_domain][i_part].n_face > 0 &&
         part_ext->parts[i_domain][i_part].face_ln_to_gn != NULL) {
        have_face = 1;
      }

    }
  }


  // En gros noeud centré avec toutes les connectivités
  if(is_describe_vtx == 1 &&
     (is_describe_edge == 0 || is_describe_face == 0) &&
     (have_edge       == 1 && have_face == 1)) {

    // Rebuild domaine_interface in distributed frame
    // PDM_domain_interface* dintrf = PDM_part_domain_interface_to_domain_interface()


    int          **pn_vtx         = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
    int          **pn_edge        = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
    int          **pn_face        = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
    PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));
    PDM_g_num_t ***pedge_ln_to_gn = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));
    PDM_g_num_t ***pface_ln_to_gn = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));
    int         ***pedge_vtx_idx  = (int         ***) malloc( part_ext->n_domain * sizeof(int         **));
    int         ***pedge_vtx      = (int         ***) malloc( part_ext->n_domain * sizeof(int         **));
    int         ***pface_edge_idx = (int         ***) malloc( part_ext->n_domain * sizeof(int         **));
    int         ***pface_edge     = (int         ***) malloc( part_ext->n_domain * sizeof(int         **));

    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      pn_vtx        [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
      pn_edge       [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
      pn_face       [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
      pvtx_ln_to_gn [i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
      pedge_ln_to_gn[i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
      pface_ln_to_gn[i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
      pedge_vtx_idx [i_domain] = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
      pedge_vtx     [i_domain] = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
      pface_edge_idx[i_domain] = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
      pface_edge    [i_domain] = (int         **) malloc( part_ext->n_domain * sizeof(int         *));

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_vtx        [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
        pn_edge       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
        pn_face       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
        pvtx_ln_to_gn [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
        pedge_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
        pface_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
        pface_edge_idx[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_edge_idx;
        pface_edge    [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_edge;
        // pedge_vtx     [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_vtx;
        pedge_vtx     [i_domain][i_part] = (int *) malloc((2 * pn_edge[i_domain][i_part]) * sizeof(int));


        int *_edge_vtx = part_ext->parts[i_domain][i_part].edge_vtx;
        for(int i_edge = 0; i_edge < pn_edge[i_domain][i_part]; ++i_edge) {
          pedge_vtx     [i_domain][i_part][2*i_edge  ] =  _edge_vtx[2*i_edge  ];
          pedge_vtx     [i_domain][i_part][2*i_edge+1] = -_edge_vtx[2*i_edge+1];
          // printf("i_edge = %i (%i)- i_vtx1 = %i | i_vtx2 = %i \n", i_edge, (int) pedge_ln_to_gn[i_domain][i_part][i_edge], _edge_vtx[2*i_edge  ], -_edge_vtx[2*i_edge+1]);
        }

        int _nedge = pn_edge[i_domain][i_part];
        pedge_vtx_idx [i_domain][i_part] = malloc((_nedge+1) * sizeof(int));

        pedge_vtx_idx [i_domain][i_part][0] = 0;
        for(int i = 0; i < _nedge; ++i) {
          pedge_vtx_idx [i_domain][i_part][i+1] = pedge_vtx_idx [i_domain][i_part][i] + 2;
        }
      }
    }

    if(is_describe_edge == 0) {

      // Translate
      printf("Translate vtx to edge ... \n");
      PDM_part_domain_interface_add(part_ext->pdi,
                                    PDM_BOUND_TYPE_VTX,
                                    PDM_BOUND_TYPE_EDGE,
                                    part_ext->n_part,
                                    pn_vtx,
                                    pvtx_ln_to_gn,
                                    pn_edge,
                                    pedge_ln_to_gn,
                                    pedge_vtx_idx,
                                    pedge_vtx,
                                    1); // Connectivity_is_signed
      printf("Translate vtx to edge END \n");
    }


    if(have_face == 1 &&  is_describe_face == 0) {

      // Translate
      printf("Translate edge to face ... \n");
      // Translate
      PDM_part_domain_interface_add(part_ext->pdi,
                                    PDM_BOUND_TYPE_EDGE,
                                    PDM_BOUND_TYPE_FACE,
                                    part_ext->n_part,
                                    pn_edge,
                                    pedge_ln_to_gn,
                                    pn_face,
                                    pface_ln_to_gn,
                                    pface_edge_idx,
                                    pface_edge,
                                    1);// Connectivity_is_signed
    }


    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(pedge_vtx_idx [i_domain][i_part]);
        free(pedge_vtx     [i_domain][i_part]);
      }
      free(pn_vtx        [i_domain]);
      free(pn_edge       [i_domain]);
      free(pn_face       [i_domain]);
      free(pvtx_ln_to_gn [i_domain]);
      free(pedge_ln_to_gn[i_domain]);
      free(pface_ln_to_gn[i_domain]);
      free(pedge_vtx_idx [i_domain]);
      free(pedge_vtx     [i_domain]);
      free(pface_edge_idx[i_domain]);
      free(pface_edge    [i_domain]);
    }
    free(pn_vtx        );
    free(pn_edge       );
    free(pn_face       );
    free(pvtx_ln_to_gn );
    free(pedge_ln_to_gn);
    free(pface_ln_to_gn);
    free(pedge_vtx_idx );
    free(pedge_vtx     );
    free(pface_edge_idx );
    free(pface_edge     );

  } else if (is_describe_face == 1) {

    // assert(is_describe_vtx == 0);

    // Faire la méthode de Julien face2vtx

  }


}


static
void
_warm_up_domain_interface
(
  PDM_part_extension_t *part_ext,
  PDM_bound_type_t      interface_kind
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  // int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    // n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  PDM_g_num_t **opp_interface_and_gnum = NULL;
  int         **current_lentity        = NULL;
  int         **current_sens           = NULL;
  int          *n_current_lentity      = NULL;

  int          **gentity2_entity1_n    = NULL;
  PDM_g_num_t  **gentity2_entity1      = NULL;

  PDM_g_num_t **entity1_opp_gnum_and_interface = NULL;
  int         **entity1_current_lentity        = NULL;
  int         **entity1_current_sens           = NULL;
  int          *n_cur_interface_entity1        = NULL;

  if(interface_kind == PDM_BOUND_TYPE_VTX) {
    assert(part_ext->opp_interface_and_gnum_vtx == NULL);
    assert(part_ext->cur_interface_vtx          == NULL);
    part_ext->opp_interface_and_gnum_vtx = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *) );
    part_ext->cur_interface_vtx          = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
    part_ext->cur_sens_vtx               = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
    part_ext->n_cur_interface_vtx        = (int          *) malloc( n_part_loc_all_domain * sizeof(int          ) );

    opp_interface_and_gnum = part_ext->opp_interface_and_gnum_vtx;
    current_lentity        = part_ext->cur_interface_vtx;
    current_sens           = part_ext->cur_sens_vtx;
    n_current_lentity      = part_ext->n_cur_interface_vtx;
  } else if( interface_kind == PDM_BOUND_TYPE_EDGE) {
    assert(part_ext->opp_interface_and_gnum_edge == NULL);
    assert(part_ext->cur_interface_edge          == NULL);
    part_ext->opp_interface_and_gnum_edge = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *) );
    part_ext->cur_interface_edge          = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
    part_ext->cur_sens_edge               = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
    part_ext->n_cur_interface_edge        = (int          *) malloc( n_part_loc_all_domain * sizeof(int          ) );


    gentity2_entity1_n = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
    gentity2_entity1   = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *) );

    opp_interface_and_gnum = part_ext->opp_interface_and_gnum_edge;
    current_lentity        = part_ext->cur_interface_edge;
    current_sens           = part_ext->cur_sens_edge;
    n_current_lentity      = part_ext->n_cur_interface_edge;
  } else if( interface_kind == PDM_BOUND_TYPE_FACE) {
    assert(part_ext->opp_interface_and_gnum_face == NULL);
    assert(part_ext->cur_interface_face          == NULL);
    part_ext->opp_interface_and_gnum_face = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *) );
    part_ext->cur_interface_face          = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
    part_ext->cur_sens_face               = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
    part_ext->n_cur_interface_face        = (int          *) malloc( n_part_loc_all_domain * sizeof(int          ) );

    gentity2_entity1_n = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
    gentity2_entity1   = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *) );

    opp_interface_and_gnum = part_ext->opp_interface_and_gnum_face;
    current_lentity        = part_ext->cur_interface_face;
    current_sens           = part_ext->cur_sens_face;
    n_current_lentity      = part_ext->n_cur_interface_face;
  }

  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
    opp_interface_and_gnum[i_part] = NULL;
    current_lentity       [i_part] = NULL;
    current_sens          [i_part] = NULL;
    n_current_lentity     [i_part] = 0;
  }

  if(part_ext->pdi != NULL) {
    int is_describe = PDM_part_domain_interface_exist_get(part_ext->pdi,
                                                         interface_kind);
    int is_describe_l = is_describe;
    PDM_MPI_Allreduce(&is_describe_l, &is_describe, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

    if(is_describe == 0) {
      if(gentity2_entity1_n != NULL) {
        free(gentity2_entity1_n);
        free(gentity2_entity1);
      }
      return;
    }
  } else {
    if(gentity2_entity1_n != NULL) {
      free(gentity2_entity1_n);
      free(gentity2_entity1);
    }
    return;
  }


  int         **neighbor_interface = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
  int         **neighbor_idx       = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
  int         **neighbor_desc      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
  int          *n_entity           = (int          *) malloc( n_part_loc_all_domain * sizeof(int          ) );
  PDM_g_num_t **entity_ln_to_gn    = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *) );
  PDM_g_num_t **entity1_ln_to_gn   = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *) );

  int **pn_entity = malloc(part_ext->n_domain * sizeof(int *));

  int connectivity_is_signed = 0;
  int shift_part   = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    pn_entity[i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int));

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      n_entity[i_part+shift_part] = 0;

      if(interface_kind == PDM_BOUND_TYPE_VTX) {
        pn_entity      [i_domain][i_part]  = part_ext->parts[i_domain][i_part].n_vtx;
        n_entity       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
        entity_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      } else if( interface_kind == PDM_BOUND_TYPE_EDGE) {
        pn_entity       [i_domain][i_part]  = part_ext->parts[i_domain][i_part].n_edge;
        n_entity        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_edge;
        entity_ln_to_gn [i_part+shift_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;

        gentity2_entity1_n[i_part+shift_part] = malloc( (pn_entity[i_domain][i_part]) * sizeof(int));
        int         *_edge_vtx    = part_ext->parts[i_domain][i_part].edge_vtx;
        PDM_g_num_t *vtx_ln_to_gn = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
        gentity2_entity1  [i_part+shift_part] = malloc( (2 * pn_entity[i_domain][i_part]) * sizeof(PDM_g_num_t));


        // Reform complete connectivity and replace implicit sens of edge_vtx
        for(int i = 0; i < pn_entity      [i_domain][i_part]; ++i) {
          gentity2_entity1_n[i_part+shift_part][i] = 2;
          int i_vtx1 = _edge_vtx[2*i  ];
          int i_vtx2 = _edge_vtx[2*i+1];

          gentity2_entity1  [i_part+shift_part][2*i  ] =  vtx_ln_to_gn[i_vtx1-1];
          gentity2_entity1  [i_part+shift_part][2*i+1] = -vtx_ln_to_gn[i_vtx2-1];
          // gentity2_entity1  [i_part+shift_part][2*i+1] =  vtx_ln_to_gn[i_vtx2-1]; // if connectivity_is_signed == 0
        }
        connectivity_is_signed = 1;
        // connectivity_is_signed = 0;

        entity1_opp_gnum_and_interface      = part_ext->opp_interface_and_gnum_vtx;
        entity1_current_lentity             = part_ext->cur_interface_vtx;
        entity1_current_sens                = part_ext->cur_sens_vtx;
        n_cur_interface_entity1             = part_ext->n_cur_interface_vtx;
        entity1_ln_to_gn[i_part+shift_part] = vtx_ln_to_gn;

      } else if( interface_kind == PDM_BOUND_TYPE_FACE) {
        pn_entity      [i_domain][i_part]  = part_ext->parts[i_domain][i_part].n_face;
        n_entity       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
        entity_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;

        int          _n_face        = part_ext->parts[i_domain][i_part].n_face;

        // 2 choices : defined from edge or direclty for vtx
        if(part_ext->parts[i_domain][i_part].n_edge != 0) {
          int         *_face_edge_idx = part_ext->parts[i_domain][i_part].face_edge_idx;
          int         *_face_edge     = part_ext->parts[i_domain][i_part].face_edge;
          PDM_g_num_t *edge_ln_to_gn  = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
          assert(_face_edge_idx != NULL);
          assert(_face_edge     != NULL);

          gentity2_entity1_n[i_part+shift_part] = malloc( _n_face * sizeof(int));
          gentity2_entity1  [i_part+shift_part] = malloc( _face_edge_idx[_n_face] * sizeof(PDM_g_num_t));

          for(int i = 0; i < _n_face; ++i) {
            gentity2_entity1_n[i_part+shift_part][i] = _face_edge_idx[i+1] - _face_edge_idx[i];
            for(int idx_edge = _face_edge_idx[i]; idx_edge < _face_edge_idx[i+1]; ++idx_edge) {
              int i_edge = PDM_ABS (_face_edge[idx_edge])-1;
              int sgn    = PDM_SIGN(_face_edge[idx_edge]);
              gentity2_entity1  [i_part+shift_part][idx_edge] = sgn * edge_ln_to_gn[i_edge];
            }
          }
          connectivity_is_signed = 1;

          entity1_opp_gnum_and_interface      = part_ext->opp_interface_and_gnum_edge;
          entity1_current_lentity             = part_ext->cur_interface_edge;
          entity1_current_sens                = part_ext->cur_sens_edge;
          n_cur_interface_entity1             = part_ext->n_cur_interface_edge;
          entity1_ln_to_gn[i_part+shift_part] = edge_ln_to_gn;

        } else {
          assert(part_ext->parts[i_domain][i_part].face_vtx != NULL);
          int         *_face_vtx_idx = part_ext->parts[i_domain][i_part].face_vtx_idx;
          int         *_face_vtx     = part_ext->parts[i_domain][i_part].face_vtx;
          PDM_g_num_t *vtx_ln_to_gn  = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;

          gentity2_entity1_n[i_part+shift_part] = malloc( _n_face * sizeof(int));
          gentity2_entity1  [i_part+shift_part] = malloc( _face_vtx_idx[_n_face] * sizeof(PDM_g_num_t));

          for(int i = 0; i < _n_face; ++i) {
            gentity2_entity1_n[i_part+shift_part][i] = _face_vtx_idx[i+1] - _face_vtx_idx[i];
            for(int idx_vtx = _face_vtx_idx[i]; idx_vtx < _face_vtx_idx[i+1]; ++idx_vtx) {
              int i_vtx = PDM_ABS (_face_vtx[idx_vtx])-1;
              gentity2_entity1  [i_part+shift_part][idx_vtx] = vtx_ln_to_gn[i_vtx];
            }
          }
          connectivity_is_signed = 0;

          entity1_opp_gnum_and_interface      = part_ext->opp_interface_and_gnum_vtx;
          entity1_current_lentity             = part_ext->cur_interface_vtx;
          entity1_current_sens                = part_ext->cur_sens_vtx;
          n_cur_interface_entity1             = part_ext->n_cur_interface_vtx;
          entity1_ln_to_gn[i_part+shift_part] = vtx_ln_to_gn;

        }
      }

    }

    shift_part   += part_ext->n_part              [i_domain];
  }

  int         **pdi_neighbor_idx         = NULL;
  int         **pdi_neighbor             = NULL;
  int           n_composed_interface     = 0;
  int          *composed_interface_idx   = NULL;
  int          *composed_interface       = NULL;
  PDM_g_num_t  *composed_ln_to_gn_sorted = NULL;

  /*
   * Attention ici car pour différentes entités on obtient des nouvelles interfaces !!!!
   * Je pense q'il faudrait finalement avoir part_ext->composed_interface + part_ext->composed_interface_ln_to_gn
   * Car sinon on a des numero absolu vtx / face par exemple qui coîncide pas / n'ont pas de sens
   */
  PDM_part_domain_interface_as_graph(part_ext->pdi,
                                     interface_kind,
                                     pn_entity,
                                     NULL,
                                     &pdi_neighbor_idx,
                                     &pdi_neighbor,
                                     &n_composed_interface,
                                     &composed_interface_idx,
                                     &composed_interface,
                                     &composed_ln_to_gn_sorted);

  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
    neighbor_idx      [i_part] = pdi_neighbor_idx[i_part];
    neighbor_desc     [i_part] = malloc( 3 * (neighbor_idx [i_part][n_entity[i_part]]) * sizeof(int));
    neighbor_interface[i_part] = malloc(     (neighbor_idx [i_part][n_entity[i_part]]) * sizeof(int));

    /* Copy */
    for(int i = 0; i < neighbor_idx [i_part][n_entity[i_part]]; ++i) {
      neighbor_desc     [i_part][3*i  ] = pdi_neighbor[i_part][4*i  ];
      neighbor_desc     [i_part][3*i+1] = pdi_neighbor[i_part][4*i+1];
      neighbor_desc     [i_part][3*i+2] = pdi_neighbor[i_part][4*i+2];
      neighbor_interface[i_part][  i  ] = pdi_neighbor[i_part][4*i+3];
    }
    // PDM_log_trace_graph_nuplet_int(neighbor_idx[i_part], neighbor_desc[i_part], 3, n_entity[i_part], "neighbor_desc (debug) :");
    // PDM_log_trace_graph_nuplet_int(neighbor_idx[i_part], pdi_neighbor[i_part], 4, n_entity[i_part], "neighbor_desc (debug) :");
    free(pdi_neighbor[i_part]);

  }
  free(pdi_neighbor_idx);
  free(pdi_neighbor);
  free(composed_interface_idx);
  free(composed_interface);
  free(composed_ln_to_gn_sorted);

  /*
   * All partition of all domain is treated in the same time
   *    -> All partition are flatten along i_domain
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           n_entity,
                                                           neighbor_idx,
                                                           neighbor_desc);

  PDM_g_num_t **entity_ln_to_gn_opp = NULL;
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                  (void **) entity_ln_to_gn,
                            NULL,
                 (void ***)&entity_ln_to_gn_opp);

  // Exchange downside connectivity to know the incoming sens of the exchange entity
  int         **entity2_entity1_opp_n = NULL;
  PDM_g_num_t **entity2_entity1_opp   = NULL;
  if(gentity2_entity1 != NULL) {
    PDM_distant_neighbor_exch(dn,
                              sizeof(PDM_g_num_t),
                              PDM_STRIDE_VAR_INTERLACED,
                              -1,
                              gentity2_entity1_n,
                    (void **) gentity2_entity1,
                             &entity2_entity1_opp_n,
                   (void ***)&entity2_entity1_opp);
  }


  shift_part   = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    // int n_part_total = part_ext->n_tot_part_by_domain[i_domain];

    /* First loop to count */
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      // int         *_neighbor_n          = neighbor_n         [i_part+shift_part];
      int         *_neighbor_idx        = neighbor_idx       [i_part+shift_part];
      int         *_neighbor_interface  = neighbor_interface [i_part+shift_part];
      PDM_g_num_t *_entity_ln_to_gn     = entity_ln_to_gn    [i_part+shift_part];
      PDM_g_num_t *_entity_ln_to_gn_opp = entity_ln_to_gn_opp[i_part+shift_part];

      opp_interface_and_gnum[i_part+shift_part] = malloc( 2 * _neighbor_idx[n_entity[i_part+shift_part]] * sizeof(PDM_g_num_t));
      current_lentity       [i_part+shift_part] = malloc(     _neighbor_idx[n_entity[i_part+shift_part]] * sizeof(int        ));
      current_sens          [i_part+shift_part] = malloc(     _neighbor_idx[n_entity[i_part+shift_part]] * sizeof(int        ));
      PDM_g_num_t *_opp_interface_and_gnum = opp_interface_and_gnum[i_part+shift_part];
      int         *_current_lentity = current_lentity[i_part+shift_part];
      int         *_current_sens    = current_sens   [i_part+shift_part];

      /* For each border count number of opposite */
      int idx_write     = 0;
      int idx_read      = 0;
      int idx_read_recv = 0;
      for(int i_entity = 0; i_entity < n_entity[i_part+shift_part]; ++i_entity) {

        /*
         * Treatement of sens
         */
        int idx_write2 = idx_write;
        if(gentity2_entity1 != NULL) {
          // log_trace("-------------------------- i_entity = %i \n", i_entity);
          for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {
            int strid = entity2_entity1_opp_n[i_part+shift_part][idx_entity];
            assert(strid == gentity2_entity1_n[i_part+shift_part][i_entity]);

            /*
             * Translate all incoming data in current partition info
             */
            // log_trace(" \t \t i_entity = %i %i \n", _entity_ln_to_gn[i_entity], _entity_ln_to_gn_opp[idx_entity]);

            // Dans le cas ou on a une entité qui rebondit sur nous on ignore car gérer par les raccord de partitionnement
            if(_entity_ln_to_gn[i_entity] == _entity_ln_to_gn_opp[idx_entity]) {
              idx_read_recv += strid;
              continue; // Cause we are the same entities, no interface need to be manage
            }

            for(int p = 0; p < strid; ++p) {
              PDM_g_num_t gopp    = PDM_ABS (entity2_entity1_opp[i_part+shift_part][idx_read_recv+p]);
              int         sgn_opp = PDM_SIGN(entity2_entity1_opp[i_part+shift_part][idx_read_recv+p]);

              PDM_g_num_t search_elmt[2] = {gopp, _neighbor_interface[idx_entity]};
              // printf("Search elmt = %i %i\n", gopp, _neighbor_interface[idx_entity]);
              // printf("n_cur_interface_entity1[i_part+shift_part] = %i \n", n_cur_interface_entity1[i_part+shift_part]);
              // PDM_log_trace_array_long(entity1_opp_gnum_and_interface[i_part+shift_part], 2 * n_cur_interface_entity1[i_part+shift_part], "entity1_opp_gnum_and_interface ::");

              // log_trace("i_entity = %i | idx_entity = %i \n", i_entity, idx_entity);
              // log_trace("Search elmt = %i %i \n", gopp, _neighbor_interface[idx_entity]);
              int pos = PDM_order_binary_search_long(search_elmt, entity1_opp_gnum_and_interface[i_part+shift_part], 2, n_cur_interface_entity1[i_part+shift_part]);
              // log_trace("pos = %i \n", pos);
              assert(pos != -1);

              int         lentity        = entity1_current_lentity[i_part+shift_part][pos    ];
              PDM_g_num_t gopp_translate = entity1_ln_to_gn       [i_part+shift_part][lentity];
              int         sgn_translate  = entity1_current_sens   [i_part+shift_part][pos    ];

              entity2_entity1_opp[i_part+shift_part][idx_read_recv+p] = sgn_translate * sgn_opp * gopp_translate;

            }

            // log_trace("current connectivity =");
            // for(int k = 0; k < strid; ++k) {
            //   log_trace("%i ", gentity2_entity1[i_part+shift_part][idx_read+k]);
            // }
            // log_trace("\n");


            // log_trace("\t --> ");
            // for(int k = 0; k < strid; ++k) {
            //   log_trace("%i ", entity2_entity1_opp[i_part+shift_part][idx_read_recv+k]);
            // }
            // log_trace("\n");

            /*
             * Determine sens
             */
            int sens = 0;
            if(connectivity_is_signed == 0) {

              if(strid == 2) { // Cas particulier edge

                PDM_g_num_t i1 = gentity2_entity1[i_part+shift_part][idx_read  ];
                PDM_g_num_t i2 = gentity2_entity1[i_part+shift_part][idx_read+1];

                PDM_g_num_t j1 = entity2_entity1_opp[i_part+shift_part][idx_read_recv  ];
                PDM_g_num_t j2 = entity2_entity1_opp[i_part+shift_part][idx_read_recv+1];

                if(i1 == j1) {
                  sens = 1;
                } else {
                  assert(i1 == j2);
                  assert(i2 == j1);
                  sens = -1;
                }
              } else {

                int idx_cur_min = 0;
                int idx_opp_min = 0;
                PDM_g_num_t gcur = gentity2_entity1   [i_part+shift_part][idx_read];
                PDM_g_num_t gopp = entity2_entity1_opp[i_part+shift_part][idx_read_recv];

                for(int k = 0; k < strid; ++k) {
                  if(gentity2_entity1[i_part+shift_part][idx_read+k] < gcur) {
                    gcur = gentity2_entity1   [i_part+shift_part][idx_read+k];
                    idx_cur_min = k;
                  }
                  if(entity2_entity1_opp[i_part+shift_part][idx_read_recv+k] < gopp) {
                    gopp = entity2_entity1_opp[i_part+shift_part][idx_read_recv+k];
                    idx_opp_min = k;
                  }
                }

                int next_idx_cur_min = (idx_cur_min + 1         ) % strid;
                int prev_idx_cur_min = (idx_cur_min - 1  + strid) % strid;
                int next_idx_opp_min = (idx_opp_min + 1         ) % strid;
                // int prev_idx_opp_min = (idx_opp_min-1) % strid;

                if(gentity2_entity1[i_part+shift_part][idx_read+next_idx_cur_min] == entity2_entity1_opp[i_part+shift_part][idx_read_recv+next_idx_opp_min])  {
                  sens = 1;
                } else {
                  // log_trace("strid      = %i\n", strid);
                  // log_trace("idx_cur_min      = %i\n", idx_cur_min);
                  // log_trace("next_idx_cur_min = %i\n", next_idx_cur_min);
                  // log_trace("prev_idx_cur_min = %i\n", prev_idx_cur_min);
                  // log_trace("next_idx_opp_min = %i\n", next_idx_opp_min);
                  // log_trace("val1 = %i| val2 = %i \n", gentity2_entity1[i_part+shift_part][idx_read+prev_idx_cur_min], entity2_entity1_opp[i_part+shift_part][idx_read_recv+next_idx_opp_min]);
                  assert(gentity2_entity1[i_part+shift_part][idx_read+prev_idx_cur_min] == entity2_entity1_opp[i_part+shift_part][idx_read_recv+next_idx_opp_min]);
                  sens = -1;
                }
              }
            } else {
              for(int k = 0; k < strid; ++k) {
                PDM_g_num_t gcur    = PDM_ABS (gentity2_entity1[i_part+shift_part][idx_read+k]);
                int         sgn_cur = PDM_SIGN(gentity2_entity1[i_part+shift_part][idx_read+k]);
                int lsens = 0;
                for(int p = 0; p < strid; ++p) {
                  PDM_g_num_t gopp   = PDM_ABS (entity2_entity1_opp[i_part+shift_part][idx_read_recv+p]);

                  if(gopp == gcur) {
                    int sgn_opp   = PDM_SIGN(entity2_entity1_opp[i_part+shift_part][idx_read_recv+p]);
                    if(sgn_opp != sgn_cur) { lsens = -1;}
                    else                   { lsens =  1;}
                  }

                  // In realease  this test is mandatory
                  if(lsens != 0) {
                    break;
                  }
                }

                assert(lsens != 0);

                if(sens != 0) {
                  assert(lsens == sens);
                }
                sens = lsens;

              }
            }
            assert(sens != 0);


            if(_entity_ln_to_gn_opp[idx_entity] != _entity_ln_to_gn[i_entity]) {
              _current_sens   [idx_write2++] = sens;
            }

            idx_read_recv += strid;
          }
          idx_read += gentity2_entity1_n[i_part+shift_part][i_entity];
        } /* END gentity2_entity1 != NULL */


        for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {

          // log_trace("(%i) _entity_ln_to_gn_opp[%i] = %i \n", _entity_ln_to_gn[i_entity], idx_entity, _entity_ln_to_gn_opp[idx_entity]);
          // _opp_interface_and_gnum[2*idx_entity  ] =  _entity_ln_to_gn_opp[idx_entity];
          // _opp_interface_and_gnum[2*idx_entity+1] =  _neighbor_interface [idx_entity]; // Opposite interface so revert sign
          // _current_lentity[idx_entity] = i_entity;

          if(_entity_ln_to_gn_opp[idx_entity] != _entity_ln_to_gn[i_entity]) {
            _opp_interface_and_gnum[2*idx_write  ] = _entity_ln_to_gn_opp[idx_entity];
            _opp_interface_and_gnum[2*idx_write+1] = _neighbor_interface[idx_entity];
            _current_lentity[idx_write] = i_entity;
            if(gentity2_entity1 == NULL) {
              _current_sens[idx_write] = 1;
            }
            idx_write++;
          }
        }
      }

      /*
       * Tri indirect sur le doublet et on order _current_lentity
       */
      // idx_write = _neighbor_idx[n_entity[i_part+shift_part]];
      int* order = malloc(     _neighbor_idx[n_entity[i_part+shift_part]] * sizeof(int        ));
      PDM_order_gnum_s(_opp_interface_and_gnum, 2, order, idx_write);

      PDM_g_num_t *unique_opp_interface_and_gnum = malloc( 2 * _neighbor_idx[n_entity[i_part+shift_part]] * sizeof(PDM_g_num_t));
      int         *unique_current_lentity        = malloc(     _neighbor_idx[n_entity[i_part+shift_part]] * sizeof(int        ));
      int         *unique_current_sens           = malloc(     _neighbor_idx[n_entity[i_part+shift_part]] * sizeof(int        ));
      // PDM_log_trace_array_long(opp_interface_and_gnum[i_part+shift_part], 2  * idx_write , "_opp_interface_and_gnum (Avant unique) ");

      int n_unique = 0;

      PDM_g_num_t last_elmt = -1;
      PDM_g_num_t last_inte = -4000000;
      for(int i = 0; i < idx_write; ++i) {
        int old_order = order[i];
        PDM_g_num_t curr_elmt   = _opp_interface_and_gnum[2*old_order  ];
        PDM_g_num_t curr_inte   = _opp_interface_and_gnum[2*old_order+1];

        int is_same = _is_same_doublet(last_elmt, last_inte, curr_elmt, curr_inte);
        if(is_same == 0) {
          unique_opp_interface_and_gnum[2*n_unique  ] = curr_elmt;
          unique_opp_interface_and_gnum[2*n_unique+1] = curr_inte;
          unique_current_lentity[n_unique] = _current_lentity[old_order];
          unique_current_sens   [n_unique] = _current_sens   [old_order];
          n_unique++;
          last_elmt = curr_elmt;
          last_inte = curr_inte;
        }

      }

      free(opp_interface_and_gnum[i_part+shift_part]);
      free(current_lentity       [i_part+shift_part]);
      free(current_sens          [i_part+shift_part]);

      unique_opp_interface_and_gnum = realloc(unique_opp_interface_and_gnum, 2 * n_unique * sizeof(PDM_g_num_t));
      unique_current_lentity        = realloc(unique_current_lentity       ,     n_unique * sizeof(int        ));
      unique_current_sens           = realloc(unique_current_sens          ,     n_unique * sizeof(int        ));

      opp_interface_and_gnum[i_part+shift_part] = unique_opp_interface_and_gnum;
      current_lentity       [i_part+shift_part] = unique_current_lentity;
      current_sens          [i_part+shift_part] = unique_current_sens   ;

      free(order);

      n_current_lentity[i_part+shift_part] = n_unique; //_neighbor_idx[n_entity[i_part+shift_part]];

      if(0 == 1) {
        PDM_log_trace_array_long(opp_interface_and_gnum[i_part+shift_part], 2  * n_current_lentity[i_part+shift_part] , "_opp_interface_and_gnum : ");
        PDM_log_trace_array_int(current_lentity        [i_part+shift_part],      n_current_lentity[i_part+shift_part] , "_current_lentity        : ");
        PDM_log_trace_array_int(current_sens           [i_part+shift_part],      n_current_lentity[i_part+shift_part] , "_current_sens           : ");
      }

      free(entity_ln_to_gn_opp[i_part+shift_part]);


    }
    shift_part   += part_ext->n_part              [i_domain];
  }

  PDM_distant_neighbor_free(dn);
  free(entity_ln_to_gn_opp);

  if(gentity2_entity1 != NULL) {
    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
      free(gentity2_entity1_n    [i_part]);
      free(gentity2_entity1      [i_part]);
      free(entity2_entity1_opp_n[i_part]);
      free(entity2_entity1_opp  [i_part]);
    }
    free(gentity2_entity1_n    );
    free(gentity2_entity1      );
    free(entity2_entity1_opp_n);
    free(entity2_entity1_opp  );
  }
  free(entity1_ln_to_gn);




  /*
   * Free
   */
  shift_part   = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    /* First loop to count */
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(neighbor_idx      [i_part+shift_part]);
      free(neighbor_desc     [i_part+shift_part]);
      free(neighbor_interface[i_part+shift_part]);
    }
    free(pn_entity[i_domain]);

    shift_part   += part_ext->n_part              [i_domain];
  }
  free(pn_entity);
  free(neighbor_idx );
  free(neighbor_desc);
  free(neighbor_interface);
  free(entity_ln_to_gn);
  free(n_entity);

}



static
void
_create_cell_graph_comm
(
  PDM_part_extension_t *part_ext
)
{
  /*
   * The idea is to use the first communicator graph to exchange directly the entity_cell
   *     -> We deduced first from this graph the data necessary for distant_neigbor
   *     -> We apply exchange with the connectivity entitiy_cell
   *        The information in graph comm is usefull also for directly adressing the results
   */
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  // assert(part_ext->n_domain == 1);

  /* Si multidomain on fait un shift et tt roule */
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }
  // printf(" n_part_loc_all_domain = %i \n", n_part_loc_all_domain);

  /* We flat all partition */
  assert(part_ext->neighbor_idx   == NULL);
  assert(part_ext->neighbor_idx   == NULL);
  assert(part_ext->n_entity_bound == NULL);

  part_ext->neighbor_idx       = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->neighbor_desc      = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->neighbor_interface = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->n_entity_bound     = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ) );

  part_ext->entity_cell_opp_idx = (int  **) malloc( n_part_loc_all_domain * sizeof(int  *) );

  part_ext->dist_neighbor_cell_n              = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->dist_neighbor_cell_idx            = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->dist_neighbor_cell_desc           = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->dist_neighbor_cell_interface      = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->unique_order_dist_neighbor_cell   = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  part_ext->n_unique_order_dist_neighbor_cell = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ) );


  // Get connectivity between domain as a graph
  int **pdi_neighbor_idx = NULL;
  int **pdi_neighbor     = NULL;

  assert(part_ext->n_composed_interface     == 0);
  assert(part_ext->composed_interface_idx   == NULL);
  assert(part_ext->composed_interface       == NULL);
  assert(part_ext->composed_ln_to_gn_sorted == NULL);

  if(part_ext->pdi != NULL) {
    if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
      int **pn_face = malloc(part_ext->n_domain * sizeof(int *));

      for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
        pn_face[i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int));
        for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
          pn_face[i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
        }
      }

      PDM_part_domain_interface_as_graph(part_ext->pdi,
                                         PDM_BOUND_TYPE_FACE,
                                         pn_face,
                                         NULL,
                                         &pdi_neighbor_idx,
                                         &pdi_neighbor,
                                         &part_ext->n_composed_interface,
                                         &part_ext->composed_interface_idx,
                                         &part_ext->composed_interface,
                                         &part_ext->composed_ln_to_gn_sorted);

      for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
        free(pn_face[i_domain]);
      }
      free(pn_face);

    } else if (part_ext->extend_type == PDM_EXTEND_FROM_VTX){

      int **pn_vtx = malloc(part_ext->n_domain * sizeof(int *));
      for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
        pn_vtx[i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int));
        for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
          pn_vtx[i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
        }
      }

      PDM_part_domain_interface_as_graph(part_ext->pdi,
                                         PDM_BOUND_TYPE_VTX,
                                         pn_vtx,
                                         NULL,
                                         &pdi_neighbor_idx,
                                         &pdi_neighbor,
                                         &part_ext->n_composed_interface,
                                         &part_ext->composed_interface_idx,
                                         &part_ext->composed_interface,
                                         &part_ext->composed_ln_to_gn_sorted);

      for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
        free(pn_vtx[i_domain]);
      }
      free(pn_vtx);
    } else {
      PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_compute wrong extend_type \n");
    }

  }

  // Begin with exchange by the connectivity the cell opposite number
  int shift_part   = 0;
  int shift_part_g = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    int n_part_total = part_ext->n_tot_part_by_domain[i_domain];

    /* First loop to count */
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_part_entity_bound_tot = -1;

      // printf(" i_part+shift_part = %i \n", i_part+shift_part);
      int* _entity_part_bound = NULL;
      if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        // part_ext->n_entity_bound[i_part+shift_part] = n_part_entity_bound_tot;
        n_part_entity_bound_tot = part_ext->parts[i_domain][i_part].face_part_bound_part_idx[n_part_total];
        part_ext->n_entity_bound[i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
        _entity_part_bound = part_ext->parts[i_domain][i_part].face_part_bound;
      } else if (part_ext->extend_type == PDM_EXTEND_FROM_VTX){
        n_part_entity_bound_tot = part_ext->parts[i_domain][i_part].vtx_part_bound_part_idx[n_part_total];
        part_ext->n_entity_bound[i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
        _entity_part_bound = part_ext->parts[i_domain][i_part].vtx_part_bound;
      } else {
        PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_compute wrong extend_type \n");
      }

      part_ext->neighbor_idx       [i_part+shift_part] = (int *) malloc(    (part_ext->n_entity_bound[i_part+shift_part]+1) * sizeof(int) );

      int n_cell = part_ext->parts[i_domain][i_part].n_cell;
      part_ext->n_cell[i_part+shift_part] = n_cell;

      int* _neighbor_idx  = part_ext->neighbor_idx [i_part+shift_part];

      int* _neighbor_n = PDM_array_zeros_int(part_ext->n_entity_bound[i_part+shift_part]);

      /* Just copy to remove the 4 indices to 3 indices */
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity = _entity_part_bound[4*idx_entity]-1;
        _neighbor_n[i_entity] += 1;
      }

      /* Join between domain */
      // PDM_bound_type_t interface_kind = PDM_BOUND_TYPE_MAX;
      // if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
      //   interface_kind = PDM_BOUND_TYPE_FACE;
      // } else if (part_ext->extend_type == PDM_EXTEND_FROM_VTX){
      //   interface_kind = PDM_BOUND_TYPE_VTX;
      // } else {
      //   PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_compute wrong extend_type \n");
      // }

      // int            n_interface        = 0;
      // int           *interface_pn       = NULL;
      // PDM_g_num_t  **interface_ln_to_gn = NULL;
      // int          **interface_sgn      = NULL;
      // int          **interface_ids      = NULL;
      // int          **interface_ids_idx  = NULL;
      // int          **interface_dom      = NULL;
      // if(part_ext->pdi != NULL) {
      //   PDM_part_domain_interface_get(part_ext->pdi,
      //                                 interface_kind,
      //                                 i_domain,
      //                                 i_part,
      //                                 &interface_pn,
      //                                 &interface_ln_to_gn,
      //                                 &interface_sgn,
      //                                 &interface_ids,
      //                                 &interface_ids_idx,
      //                                 &interface_dom);
      //   n_interface = PDM_part_domain_interface_n_interface_get(part_ext->pdi);
      // }

      // /*
      //  * First step : Count interface to add in distant neighbor due to connectivity betwenn domain
      //  */
      // for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

      //   // log_trace("-------------------------------- i_interface = %i  -------------------------------- \n", i_interface);
      //   // PDM_log_trace_array_int(interface_sgn[i_interface], interface_pn[i_interface], "interface_sgn :: ");

      //   for(int idx_entity = 0; idx_entity < interface_pn[i_interface]; ++idx_entity) {

      //     // Search the first in list that is in current part/proc
      //     // int i_proc_cur   = -1;
      //     // int i_part_cur   = -1;
      //     int i_entity_cur = -1;
      //     int found        = 0;
      //     int idx_current  = -1;
      //     for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
      //       int i_proc_opp   = interface_ids[i_interface][3*j  ];
      //       int i_part_opp   = interface_ids[i_interface][3*j+1];
      //       int i_entity_opp = interface_ids[i_interface][3*j+2];

      //       if(i_proc_opp == i_rank && i_part_opp == i_part) {
      //         // i_proc_cur   = i_proc_opp;
      //         // i_part_cur   = i_part_opp;
      //         i_entity_cur = i_entity_opp;
      //         idx_current  = j;
      //         assert(found == 0);
      //         found = 1;
      //         break;
      //       }
      //     }

      //     if(!found) {
      //       continue;
      //     }

      //     // Il manque une notion de direction sinon on sait pas dans quelle sens va le raccord

      //     assert(found == 1);

      //     // log_trace("i_proc_cur = %i | i_part_cur = %i | i_entity_cur = %i | sgn = %i \n", i_proc_cur, i_part_cur, i_entity_cur, interface_sgn[i_interface][idx_entity]);

      //     // Only add the opposite part of the graph
      //     for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
      //       // int i_proc_opp   = interface_ids[i_interface][3*j  ];
      //       // int i_part_opp   = interface_ids[i_interface][3*j+1];
      //       // int i_entity_opp = interface_ids[i_interface][3*j+2];

      //       if(idx_current != j) {
      //         // log_trace("\t i_proc_opp = %i | i_part_opp = %i | i_entity_opp = %i \n", i_proc_opp, i_part_opp, i_entity_opp);
      //         _neighbor_n[i_entity_cur] += 1;
      //       }
      //     }
      //   }
      // }

      /*
       * Count comming from interface
       */
      if(pdi_neighbor_idx != NULL) {
        for(int i_entity = 0; i_entity < part_ext->n_entity_bound[i_part+shift_part]; ++i_entity) {
          _neighbor_n[i_entity] += pdi_neighbor_idx[i_part+shift_part][i_entity+1] - pdi_neighbor_idx[i_part+shift_part][i_entity];
        }
      }


      /* Compute index */
      _neighbor_idx[0] = 0;
      for(int i_entity = 0; i_entity < part_ext->n_entity_bound[i_part+shift_part]; ++i_entity) {
        _neighbor_idx[i_entity+1] = _neighbor_idx[i_entity] + _neighbor_n[i_entity];
        _neighbor_n[i_entity] = 0;
      }


      part_ext->neighbor_desc     [i_part+shift_part] = (int *) malloc( 3 * _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]] * sizeof(int) );
      part_ext->neighbor_interface[i_part+shift_part] = (int *) malloc(     _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]] * sizeof(int) );
      int* _neighbor_desc      = part_ext->neighbor_desc     [i_part+shift_part];
      int* _neighbor_interface = part_ext->neighbor_interface[i_part+shift_part];

      /* Second loop for fill */
      /* Just copy to remove the 4 indices to 3 indices */
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity = _entity_part_bound[4*idx_entity]-1;
        int idx_write = _neighbor_idx[i_entity] + _neighbor_n[i_entity]++;
        _neighbor_desc[3*idx_write  ] = _entity_part_bound[4*idx_entity+1];                // i_proc_opp;
        _neighbor_desc[3*idx_write+1] = _entity_part_bound[4*idx_entity+2]+shift_part_g-1; // i_part_opp
        _neighbor_desc[3*idx_write+2] = _entity_part_bound[4*idx_entity+3]-1;              // i_entity_opp
        _neighbor_interface[idx_write] = -40000;
      }

      /*
       * Add graph du to connectivy between domain
       */
      // for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      //   for(int idx_entity = 0; idx_entity < interface_pn[i_interface]; ++idx_entity) {

      //     // Search the first in list that is in current part/proc
      //     // int i_proc_cur   = -1;
      //     // int i_part_cur   = -1;
      //     int i_entity_cur = -1;
      //     int found        = 0;
      //     int idx_current  = -1;
      //     for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
      //       int i_proc_opp   = interface_ids[i_interface][3*j  ];
      //       int i_part_opp   = interface_ids[i_interface][3*j+1];
      //       int i_entity_opp = interface_ids[i_interface][3*j+2];

      //       if(i_proc_opp == i_rank && i_part_opp == i_part) {
      //         // i_proc_cur   = i_proc_opp;
      //         // i_part_cur   = i_part_opp;
      //         i_entity_cur = i_entity_opp;
      //         idx_current  = j;
      //         assert(found == 0);
      //         found = 1;
      //         break;
      //       }
      //     }

      //     if(!found) {
      //       continue;
      //     }

      //     assert(found == 1);

      //     // Only add the opposite part of the graph
      //     for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
      //       int i_proc_opp   = interface_ids[i_interface][3*j  ];
      //       int i_part_opp   = interface_ids[i_interface][3*j+1];
      //       int i_entity_opp = interface_ids[i_interface][3*j+2];

      //       if(idx_current != j) {
      //         // log_trace("\t i_proc_opp = %i | i_part_opp = %i | i_entity_opp = %i \n", i_proc_opp, i_part_opp, i_entity_opp);
      //         int idx_write = _neighbor_idx[i_entity_cur] + _neighbor_n[i_entity_cur]++;
      //         _neighbor_desc[3*idx_write  ] = i_proc_opp;              // i_proc_opp;
      //         _neighbor_desc[3*idx_write+1] = i_part_opp+shift_part_g; // i_part_opp
      //         _neighbor_desc[3*idx_write+2] = i_entity_opp;            // i_entity_opp
      //         _neighbor_interface[idx_write] = (i_interface+1) * interface_sgn[i_interface][idx_entity];
      //         // _neighbor_interface[idx_write] = (i_interface+1);
      //       }
      //     }
      //   }
      // }
      if(pdi_neighbor_idx != NULL) {
        for(int i_entity = 0; i_entity < part_ext->n_entity_bound[i_part+shift_part]; ++i_entity) {
          for(int idx = pdi_neighbor_idx[i_part+shift_part][i_entity]; idx < pdi_neighbor_idx[i_part+shift_part][i_entity+1]; ++idx) {
            int idx_write = _neighbor_idx[i_entity] + _neighbor_n[i_entity]++;
            _neighbor_desc[3*idx_write  ]  = pdi_neighbor[i_part+shift_part][4*idx  ];
            _neighbor_desc[3*idx_write+1]  = pdi_neighbor[i_part+shift_part][4*idx+1];
            _neighbor_desc[3*idx_write+2]  = pdi_neighbor[i_part+shift_part][4*idx+2];
            _neighbor_interface[idx_write] = pdi_neighbor[i_part+shift_part][4*idx+3];
          }
        }
      }

      free(_neighbor_n);

      if(0 == 1) {
        PDM_log_trace_array_int(_neighbor_idx , part_ext->n_entity_bound[i_part+shift_part]+1, "_neighbor_idx::");
        PDM_log_trace_array_int(_neighbor_desc , 3 * _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]], "_neighbor_desc::");
        PDM_log_trace_array_int(_neighbor_interface , _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]], "_neighbor_interface::");
      }


    }

    shift_part   += part_ext->n_part              [i_domain];
    shift_part_g += part_ext->n_tot_part_by_domain[i_domain];
  }


  /*
   * All partition of all domain is treated in the same time
   *    -> All partition are flatten along i_domain
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           part_ext->n_entity_bound,
                                                           part_ext->neighbor_idx,
                                                           part_ext->neighbor_desc);


  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            part_ext->entity_cell_n,
                  (void **) part_ext->entity_cell,
                           &part_ext->entity_cell_opp_n,
                 (void ***)&part_ext->entity_cell_opp);
  // PDM_distant_neighbor_exch_int(dn,
  //                           sizeof(int),
  //                           PDM_STRIDE_VAR,
  //                           -1,
  //                           part_ext->entity_cell_n,
  //                 (void **) part_ext->entity_cell,
  //                          &part_ext->entity_cell_opp_n,
  //                (void ***)&part_ext->entity_cell_opp);

  /* Compute idx */
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell = part_ext->parts[i_domain][i_part].n_cell;
      int* _neighbor_idx        = part_ext->neighbor_idx       [i_part+shift_part];
      int* _neighbor_desc       = part_ext->neighbor_desc      [i_part+shift_part];
      int* _neighbor_interface  = part_ext->neighbor_interface [i_part+shift_part];

      int* _entity_cell_opp_n   = part_ext->entity_cell_opp_n  [i_part+shift_part];
      int* _entity_cell_opp     = part_ext->entity_cell_opp    [i_part+shift_part];

      int* _entity_cell_idx = part_ext->entity_cell_idx[i_part+shift_part];
      int* _entity_cell     = part_ext->entity_cell    [i_part+shift_part];

      int n_entity = part_ext->n_entity_bound[i_part+shift_part]; // n_face / n_edge / n_vtx
      int n_elmt   = _neighbor_idx[n_entity];

      // int n_cell = part_ext->parts[i_domain][i_part].n_cell;
      part_ext->entity_cell_opp_idx   [i_part+shift_part] = PDM_array_new_idx_from_sizes_int(_entity_cell_opp_n, n_elmt);
      part_ext->dist_neighbor_cell_n  [i_part+shift_part] = PDM_array_zeros_int(n_cell);
      part_ext->dist_neighbor_cell_idx[i_part+shift_part] = (int *) malloc( (n_cell+1) * sizeof(int) );

      int* _entity_cell_opp_idx    = part_ext->entity_cell_opp_idx   [i_part+shift_part];
      int* _dist_neighbor_cell_n   = part_ext->dist_neighbor_cell_n  [i_part+shift_part];
      int* _dist_neighbor_cell_idx = part_ext->dist_neighbor_cell_idx[i_part+shift_part];


      // part_ext->border_cell_list[i_part+shift_part] = malloc( n_elmt * sizeof(int));
      part_ext->border_cell_list[i_part+shift_part] = malloc( n_cell * sizeof(int));
      // for(int i = 0; i < n_elmt; ++i) {
      //   part_ext->border_cell_list[i_part+shift_part][i] = -1000;
      // }

      /* For each border count number of opposite */
      for(int i_entity = 0; i_entity < n_entity; ++i_entity) {
        for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {
          for(int idx_cell = _entity_cell_idx[i_entity]; idx_cell < _entity_cell_idx[i_entity+1]; ++idx_cell) {
            int i_cell = PDM_ABS(_entity_cell[idx_cell])-1;
            _dist_neighbor_cell_n[i_cell] += _entity_cell_opp_n[idx_entity];
          }
        }
      }

      // PDM_log_trace_array_int(_dist_neighbor_cell_n , n_cell, "_dist_neighbor_cell_n::");
      // PDM_log_trace_array_int(_entity_cell_opp_n , n_part_entity_bound_tot, "_entity_cell_opp_n::");
      // PDM_log_trace_array_int(_entity_cell_opp , n_part_entity_bound_tot, "_entity_cell_opp::");

      /* Ici il faut faire les raccords entre domaine ---> Count */

      /* Compute idx */
      int idx_indic = 0;
      _dist_neighbor_cell_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _dist_neighbor_cell_idx[i_cell+1] = _dist_neighbor_cell_idx[i_cell] + _dist_neighbor_cell_n[i_cell];
        if(_dist_neighbor_cell_n[i_cell] > 0){
          part_ext->border_cell_list[i_part+shift_part][idx_indic++] = i_cell;
        }
      }
      part_ext->border_cell_list[i_part+shift_part] = realloc(part_ext->border_cell_list[i_part+shift_part], idx_indic * sizeof(int));
      /* Because cell can be associated twice */
      part_ext->n_cell_border[i_part+shift_part] = idx_indic;

      // printf(" idx_indic = %i\n", idx_indic);
      // PDM_log_trace_array_int(part_ext->border_cell_list[i_part+shift_part]  , idx_indic  , "border_cell_list::");

      /* Reset */
      PDM_array_reset_int(_dist_neighbor_cell_n, n_cell, 0);

      /* Allocate */
      part_ext->dist_neighbor_cell_desc[i_part+shift_part] = (int * ) malloc( 4 * _dist_neighbor_cell_idx[n_cell] * sizeof(int));
      int* _dist_neighbor_cell_desc = part_ext->dist_neighbor_cell_desc[i_part+shift_part];

      for(int i_entity = 0; i_entity < n_entity; ++i_entity) {
        for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {
          int i_proc_opp   = _neighbor_desc[3*idx_entity  ];
          int i_part_opp   = _neighbor_desc[3*idx_entity+1];
          int i_interf     = _neighbor_interface[idx_entity];
          // int i_entity_opp = _neighbor_desc[3*idx_entity+2];
          for(int idx_cell = _entity_cell_idx[i_entity]; idx_cell < _entity_cell_idx[i_entity+1]; ++idx_cell) {
            int i_cell = PDM_ABS(_entity_cell[idx_cell])-1;

            for(int idx_cell_opp = _entity_cell_opp_idx[idx_entity]; idx_cell_opp < _entity_cell_opp_idx[idx_entity+1]; ++idx_cell_opp) {
              int idx_write = _dist_neighbor_cell_idx[i_cell] + _dist_neighbor_cell_n[i_cell]++;
              int i_opp_cell = PDM_ABS(_entity_cell_opp[idx_cell_opp])-1; // Commence à zero
              // printf("[%i] - _entity_cell_opp[%i] = %i | idx_write = %i  \n", i_part, idx_cell_opp,  _entity_cell_opp[idx_cell_opp], idx_write);
              // _dist_neighbor_cell_desc[3*idx_write  ] = i_proc_opp;
              // _dist_neighbor_cell_desc[3*idx_write+1] = i_part_opp;
              // _dist_neighbor_cell_desc[3*idx_write+2] = i_opp_cell;

              _dist_neighbor_cell_desc[4*idx_write  ] = i_proc_opp;
              _dist_neighbor_cell_desc[4*idx_write+1] = i_part_opp;
              _dist_neighbor_cell_desc[4*idx_write+2] = i_opp_cell;
              _dist_neighbor_cell_desc[4*idx_write+3] = i_interf; // PDM_ABS(i_interf);
              // printf("i_interf = %i \n", i_interf);
              // printf("[%i][%i] _dist_neighbor_cell[%i] = %i %i %i - i_entity = %i from ii_part_opp = %i and  i_entity_opp = %i \n",
              //        i_part, i_cell, idx_write, i_proc_opp, i_part_opp, i_opp_cell, i_entity, i_part_opp, i_entity_opp);
            }
          }
        }
      }

      int* _unique_dist_neighbor_cell_idx = NULL;
      int* _unique_dist_neighbor_cell_n   = NULL;
      int* _unique_dist_neighbor_cell     = NULL;
      _unique_quadruplet(n_cell,
                         _dist_neighbor_cell_idx,
                         _dist_neighbor_cell_desc,
                         &_unique_dist_neighbor_cell_idx,
                         &_unique_dist_neighbor_cell_n,
                         &_unique_dist_neighbor_cell);
      free(_dist_neighbor_cell_idx);
      free(_dist_neighbor_cell_desc);
      free(_dist_neighbor_cell_n);

      int *_unique_order_dist_neighbor_cell = NULL;
      int n_unique = _setup_unique_order_quadruplet(n_cell,
                                                    _unique_dist_neighbor_cell_idx,
                                                    _unique_dist_neighbor_cell,
                                                    &_unique_order_dist_neighbor_cell);

      int *_unique_dist_neighbor_cell_interface = NULL;
      quadruplet_to_triplet_and_array(_unique_dist_neighbor_cell_idx[n_cell],
                                      _unique_dist_neighbor_cell,
                                      &_unique_dist_neighbor_cell_interface,
                                      &part_ext->dist_neighbor_cell_desc[i_part+shift_part]);

      part_ext->dist_neighbor_cell_idx         [i_part+shift_part] = _unique_dist_neighbor_cell_idx;
      part_ext->dist_neighbor_cell_n           [i_part+shift_part] = _unique_dist_neighbor_cell_n;
      free(_unique_dist_neighbor_cell    );

      part_ext->unique_order_dist_neighbor_cell  [i_part+shift_part] = _unique_order_dist_neighbor_cell;
      part_ext->n_unique_order_dist_neighbor_cell[i_part+shift_part] = n_unique;

      /*
       * Storee interface number
       */
      part_ext->dist_neighbor_cell_interface[i_part+shift_part] = _unique_dist_neighbor_cell_interface;

    }
    shift_part += part_ext->n_part[i_domain];
  }

  PDM_distant_neighbor_free(dn);

  // for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
  //   free(pdi_neighbor_idx[i_part]);
  //   free(pdi_neighbor    [i_part]);
  // }
  // free(pdi_neighbor_idx);
  // free(pdi_neighbor    );
  part_ext->pdi_neighbor_idx = pdi_neighbor_idx;
  part_ext->pdi_neighbor     = pdi_neighbor;


  /*
   * Panic verbose
   */
  if(0 == 1) {
    shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        int* _neighbor_idx = part_ext->neighbor_idx                [i_part+shift_part];
        int n_elmt         = _neighbor_idx[part_ext->n_entity_bound[i_part+shift_part]];
        int n_data         = part_ext->entity_cell_opp_idx[i_part+shift_part][n_elmt];

        PDM_log_trace_array_int(part_ext->entity_cell_opp_n  [i_part+shift_part], n_elmt  , "entity_cell_opp_n::");
        PDM_log_trace_array_int(part_ext->entity_cell_opp_idx[i_part+shift_part], n_elmt+1, "entity_cell_opp_idx::");
        PDM_log_trace_array_int(part_ext->entity_cell_opp    [i_part+shift_part], n_data  , "entity_cell_opp::");

        int n_cell = part_ext->parts[i_domain][i_part].n_cell;
        int* _dist_neighbor_cell_idx  = part_ext->dist_neighbor_cell_idx      [i_part+shift_part];
        int* _dist_neighbor_cell_desc = part_ext->dist_neighbor_cell_desc     [i_part+shift_part];
        int* _dist_neighbor_cell_intf = part_ext->dist_neighbor_cell_interface[i_part+shift_part];

        log_trace("*** i_part %d ***\n", i_part);
        log_trace(" _dist_neighbor_cell_idx ---------------- \n");
        for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
          log_trace("[%i] i_cell -> %i -->  ", i_part, i_cell);
          for(int idx = _dist_neighbor_cell_idx[i_cell]; idx < _dist_neighbor_cell_idx[i_cell+1]; ++idx) {
            log_trace("(%i, %i, intrf = %i) ", _dist_neighbor_cell_desc[3*idx+1], _dist_neighbor_cell_desc[3*idx+2], _dist_neighbor_cell_intf[idx]);
          }
          log_trace("\n");
        }
        log_trace(" _dist_neighbor_cell_idx ---------------- END \n");

        PDM_log_trace_array_int(_dist_neighbor_cell_idx , n_cell+1                       , "_dist_neighbor_cell_idx::");
        // PDM_log_trace_array_int(_dist_neighbor_cell_desc, 3 * _dist_neighbor_cell_idx[n_cell], "_dist_neighbor_cell_desc::");
        PDM_log_trace_graph_nuplet_int(_dist_neighbor_cell_idx,
                                       _dist_neighbor_cell_desc,
                                       3,
                                       n_cell,
                                       "_dist_neighbor_cell_desc::");

        PDM_log_trace_array_int(_dist_neighbor_cell_intf,     _dist_neighbor_cell_idx[n_cell], "_dist_neighbor_cell_intf::");
      }
      shift_part += part_ext->n_part[i_domain];
    }
  }
  // exit(1);


}


static
void
_compute_dual_graph
(
  PDM_part_extension_t   *part_ext,
  PDM_distant_neighbor_t *dn,
  int                     i_depth
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  // printf("_compute_dual_graph : %i  \n", i_depth);

  assert(i_depth > 0);
  int** prev_cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth-1];
  int** prev_cell_cell_extended_n   = part_ext->cell_cell_extended_n  [i_depth-1];
  int** prev_cell_cell_extended     = part_ext->cell_cell_extended    [i_depth-1];
  int** prev_cell_cell_interface    = part_ext->cell_cell_interface   [i_depth-1];

  // int** next_cell_cell_extended_idx = NULL;
  int** next_cell_cell_extended_n   = NULL;
  int** next_cell_cell_extended     = NULL;

  /*
   * Exchange of the previous rank
   */
  PDM_distant_neighbor_exch(dn,
                            3 * sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            prev_cell_cell_extended_n,
                  (void **) prev_cell_cell_extended,
                           &next_cell_cell_extended_n,
                 (void ***)&next_cell_cell_extended);

  int** next_cell_cell_interface_n   = NULL;
  int** next_cell_cell_interface     = NULL;
  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            prev_cell_cell_extended_n,
                  (void **) prev_cell_cell_interface,
                           &next_cell_cell_interface_n,
                 (void ***)&next_cell_cell_interface);
   // Same as next_cell_cell_extended_n
  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(next_cell_cell_interface_n[i_part+shift_part]);
    }
    shift_part += part_ext->n_part[i_domain];
  }
  free(next_cell_cell_interface_n);

  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell        = part_ext->n_cell       [i_part+shift_part];
      int n_cell_border = part_ext->n_cell_border[i_part+shift_part];

      int* _border_cell_cell_extended_n   = next_cell_cell_extended_n[i_part+shift_part];
      int* _border_cell_cell_extended     = next_cell_cell_extended  [i_part+shift_part];
      int* _border_cell_cell_interface    = next_cell_cell_interface [i_part+shift_part];
      int* _border_cell_cell_extended_idx = (int * ) malloc( (n_cell_border+1) * sizeof(int) );

      if( 0 == 1) {
        printf("prev_cell_cell_extended :: --------------------- \n");
        for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
          if( prev_cell_cell_extended_idx[i_part+shift_part][i_cell+1] > prev_cell_cell_extended_idx[i_part+shift_part][i_cell]){
            printf("[%i] i_cell -> %i -->  ", i_part, i_cell);
          }
          for(int idx = prev_cell_cell_extended_idx[i_part+shift_part][i_cell]; idx < prev_cell_cell_extended_idx[i_part+shift_part][i_cell+1]; ++idx) {
            printf("(%i, %i, intrf = %i) ", prev_cell_cell_extended[i_part+shift_part][3*idx+1],
                                            prev_cell_cell_extended[i_part+shift_part][3*idx+2],
                                            prev_cell_cell_interface[i_part+shift_part][idx]);
          }
          if( prev_cell_cell_extended_idx[i_part+shift_part][i_cell+1] > prev_cell_cell_extended_idx[i_part+shift_part][i_cell]){
            printf("\n");
          }
        }
        printf("prev_cell_cell_extended :: --------------------- END \n");
      }

      // _border_cell_cell_extended_idx[0] = 0;
      // for(int i = 0; i < n_cell_border; ++i) {
      //   _border_cell_cell_extended_idx[i+1] = _border_cell_cell_extended_idx[i] + _border_cell_cell_extended_n[i];
      // }

      int lidx_write = 0;
      _border_cell_cell_extended_idx[0] = 0;
      int idx_read = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        int n_neight = 0;
        for(int idx_neight = part_ext->dist_neighbor_cell_idx [i_part+shift_part][i_cell]; idx_neight < part_ext->dist_neighbor_cell_idx[i_part+shift_part][i_cell+1]; ++idx_neight){
          n_neight += _border_cell_cell_extended_n[idx_neight];


          // printf(" Coming in interface = %i \n", part_ext->dist_neighbor_cell_interface[i_part+shift_part][idx_neight]);
          int cur_interface = part_ext->dist_neighbor_cell_interface[i_part+shift_part][idx_neight];

          for(int j = 0; j < _border_cell_cell_extended_n[idx_neight]; ++j ) {

            // Adapt opposite with the current one
            if(_border_cell_cell_interface[idx_read] == -40000) {
              _border_cell_cell_interface[idx_read] = cur_interface;
            } else { // Composition

              if(cur_interface == -_border_cell_cell_interface[idx_read]) {
                 _border_cell_cell_interface[idx_read] = -40000;
              }

              // _border_cell_cell_interface[idx_read] = cur_interface;
              // printf("_compute_dual_graph Composition a gerer \n");
              // abort();
              // log_trace("_border_cell_cell_interface[%i] = %i \n", idx_read, _border_cell_cell_interface[idx_read]);
              // int interf_border = _border_cell_cell_interface[idx_read];
              // _border_cell_cell_interface[idx_read] = cur_interface + part_ext->n_interface * interf_border;
            }


            idx_read++;
          }



        }
        // printf(" beg = %i | end = %i | n_neight = %i \n", part_ext->dist_neighbor_cell_idx [i_part+shift_part][i_cell], part_ext->dist_neighbor_cell_idx [i_part+shift_part][i_cell+1], n_neight );
        if(n_neight > 0) {
          _border_cell_cell_extended_idx[lidx_write+1] = _border_cell_cell_extended_idx[lidx_write] + n_neight;
          lidx_write++;
        }
      }

      for(int i = 0; i < n_cell_border; ++i) {
        _border_cell_cell_extended_n[i] = _border_cell_cell_extended_idx[i+1] - _border_cell_cell_extended_idx[i];
      }

      int *_border_cell_cell_extended_and_interface = NULL;
      triplet_to_quadruplet(_border_cell_cell_extended_idx[n_cell_border],
                            _border_cell_cell_extended,
                            _border_cell_cell_interface,
                            &_border_cell_cell_extended_and_interface);
      free(_border_cell_cell_extended     );
      free(_border_cell_cell_interface     );

      int* _unique_border_cell_cell_extended_n   = NULL;
      int* _unique_border_cell_cell_extended     = NULL;
      int* _unique_border_cell_cell_extended_idx = NULL;
      _unique_quadruplet(n_cell_border,
                         _border_cell_cell_extended_idx,
                         _border_cell_cell_extended_and_interface,
                         &_unique_border_cell_cell_extended_idx,
                         &_unique_border_cell_cell_extended_n,
                         &_unique_border_cell_cell_extended);
      free(_border_cell_cell_extended_n   );
      free(_border_cell_cell_extended_idx );
      free(_border_cell_cell_extended_and_interface);

      _border_cell_cell_extended_n   = _unique_border_cell_cell_extended_n  ;
      _border_cell_cell_extended_idx = _unique_border_cell_cell_extended_idx;
      _border_cell_cell_extended     = _unique_border_cell_cell_extended    ;
      next_cell_cell_extended_n[i_part+shift_part] = _border_cell_cell_extended_n;
      next_cell_cell_extended  [i_part+shift_part] = _border_cell_cell_extended;

      int *_unique_order_border_cell_cell_extended = NULL;
      int n_unique_border = _setup_unique_order_quadruplet(n_cell_border,
                                                           _unique_border_cell_cell_extended_idx,
                                                           _unique_border_cell_cell_extended,
                                                           &_unique_order_border_cell_cell_extended);

      quadruplet_to_triplet_and_array(_unique_border_cell_cell_extended_idx[n_cell_border],
                                      _unique_border_cell_cell_extended,
                                      &_border_cell_cell_interface,
                                      &_border_cell_cell_extended);
      free(_unique_border_cell_cell_extended);
      next_cell_cell_extended  [i_part+shift_part] = _border_cell_cell_extended;
      next_cell_cell_interface [i_part+shift_part] = _border_cell_cell_interface;


      if(0 == 1) {
        PDM_log_trace_array_int(_border_cell_cell_extended_idx, n_cell_border+1, "next_border_cell_cell_extended_idx::");
        PDM_log_trace_array_int(_border_cell_cell_extended_n  , n_cell_border,   "next_border_cell_cell_extended_n::");
        PDM_log_trace_array_int(_border_cell_cell_extended    , 3 * _border_cell_cell_extended_idx[n_cell_border], "next_border_cell_cell_extended::");

        printf("_border_cell_cell_extended :: --------------------- \n");
        for(int i = 0; i < n_cell_border; ++i) {
          int i_cell = part_ext->border_cell_list[i_part+shift_part][i];
          if( _border_cell_cell_extended_idx[i+1] > _border_cell_cell_extended_idx[i]){
            printf("i_cell -> %i -->  ", i_cell);
          }
          for(int idx = _border_cell_cell_extended_idx[i]; idx < _border_cell_cell_extended_idx[i+1]; ++idx) {
            printf("(%i, %i, intrf = %i) ", _border_cell_cell_extended[3*idx+1],
                                            _border_cell_cell_extended[3*idx+2],
                                            _border_cell_cell_interface[idx]);
          }
          if( _border_cell_cell_extended_idx[i+1] > _border_cell_cell_extended_idx[i]){
            printf("\n");
          }
        }
        printf("_border_cell_cell_extended :: --------------------- END \n");

      }

      /* Now we have the opposite border_cell_cell_extended
       *   -> We need to reput all together
       *   -> The current cell has to be also extendted but caution !!!
       *      the extension depend also of border
       */

      /*
       * Generate tag to know is a local cell is a border or not
       */
      int* idx_border_cell = PDM_array_const_int(n_cell, -1);

      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        idx_border_cell[i_cell] = idx_cell; /* Keep idx_cell for imediate loop */
      }

      /* Allocate current depth and shorcut */
      part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part] = (int *) malloc( (n_cell + 1 ) * sizeof(int));
      part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part] = (int *) malloc( (n_cell     ) * sizeof(int));

      int* _prev_cell_cell_extended_n   = prev_cell_cell_extended_n  [i_part+shift_part];
      int* _prev_cell_cell_extended_idx = prev_cell_cell_extended_idx[i_part+shift_part];
      int* _prev_cell_cell_extended     = prev_cell_cell_extended    [i_part+shift_part];
      int* _prev_cell_cell_interface    = prev_cell_cell_interface   [i_part+shift_part];

      int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part];
      int* _cell_cell_extended_n   = part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part];

      int* _cell_cell_idx = part_ext->cell_cell_idx[i_part+shift_part];
      int* _cell_cell     = part_ext->cell_cell    [i_part+shift_part];

      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell] = 0;
        _cell_cell_extended_n  [i_cell] = 0;
      }
      _cell_cell_extended_idx[n_cell] = 0;

      /* Weak hash table */
      int n_unique_dist_neighbor_cell = part_ext->n_unique_order_dist_neighbor_cell[i_part+shift_part];
      int n_unique_cell_cell_extended = part_ext->n_unique_order_cell_cell_extended[i_depth-1][i_part+shift_part];

      int *_unique_prev_cell_cell_extended = part_ext->unique_order_cell_cell_extended[i_depth-1][i_part+shift_part];

      // printf("n_unique_dist_neighbor_cell = %i \n", n_unique_dist_neighbor_cell);
      // printf("n_unique_cell_cell_extended = %i \n", n_unique_cell_cell_extended);

      int *tag_border_cell_cell_extended = malloc( n_unique_border             * sizeof(int));
      int *tag_dist_neighbor_cell        = malloc( n_unique_dist_neighbor_cell * sizeof(int));
      int *tag_prev_cell_cell_extended   = malloc( n_unique_cell_cell_extended * sizeof(int));
      int *tag_interior_cell             = malloc( n_cell                      * sizeof(int));

      int n_work = n_unique_border + n_unique_dist_neighbor_cell + n_unique_cell_cell_extended;
      int *icell_to_reset = malloc( n_work * sizeof(int)); // Suralloc but don't really have choice

      int *icell_border_to_reset = malloc( n_unique_border * sizeof(int)); // Suralloc but don't really have choice

      int n_work2 = _cell_cell_idx[n_cell];
      int *icell_to_reset_interior = malloc( n_work2 * sizeof(int)); // Suralloc but don't really have choice

      for(int i = 0; i < n_unique_border; ++i) {
        tag_border_cell_cell_extended[i] = 0;
      }
      for(int i = 0; i < n_unique_dist_neighbor_cell; ++i) {
        tag_dist_neighbor_cell[i] = 0;
      }
      for(int i = 0; i < n_unique_cell_cell_extended; ++i) {
        tag_prev_cell_cell_extended[i] = 0;
      }

      for(int i = 0; i < n_cell; ++i) {
        tag_interior_cell[i] = 0;
      }

      /* Let's go - First pass to count */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        int n_unique_loc = 0;
        int n_interior_unique_loc = 0;
        int n_border_unique_loc = 0;

        assert(_cell_cell_extended_n[i_cell] == 0); // All are sorted before

        /* From interior - we add the previous rank */
        _cell_cell_extended_n[i_cell] = _prev_cell_cell_extended_n[i_cell];
        for(int idx_neight = _prev_cell_cell_extended_idx[i_cell]; idx_neight < _prev_cell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int i_unique2 = _unique_prev_cell_cell_extended[idx_neight]; // Pas la cellule le tableau est unique !!!!!
          tag_prev_cell_cell_extended[i_unique2 ] = 1;
          icell_to_reset[n_unique_loc++ ] = i_unique2;
        }

        /* From border */
        // printf("_prev_cell_cell_extended_n[%i][%i] = %i\n", i_part, i_cell, _prev_cell_cell_extended_n[i_cell]);
        // printf("_border_cell_cell_extended_n[%i][%i] = %i\n", i_part, idx_cell, _border_cell_cell_extended_n[idx_cell]);
        assert(_border_cell_cell_extended_n[idx_cell] > 0);
        // _cell_cell_extended_n[i_cell] += _border_cell_cell_extended_n[idx_cell];

        for(int idx_neight = _border_cell_cell_extended_idx[idx_cell]; idx_neight < _border_cell_cell_extended_idx[idx_cell+1]; ++idx_neight) {
          int i_unique2 = _unique_order_border_cell_cell_extended[idx_neight];
          int is_treat2 = tag_border_cell_cell_extended          [i_unique2 ];
          if(is_treat2 == 0) {
            _cell_cell_extended_n[i_cell] += 1;
            tag_border_cell_cell_extended[i_unique2]  = 1;
            icell_border_to_reset[n_border_unique_loc++] = i_unique2;
          }
        }

        /* Now we have to extend the interior */
        for(int idx_neight = _prev_cell_cell_extended_idx[i_cell]; idx_neight < _prev_cell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int i_rank_neight = _prev_cell_cell_extended[3*idx_neight  ];
          int i_part_neight = _prev_cell_cell_extended[3*idx_neight+1];
          /* We add stencil only if it's local */
          if(i_part+shift_part == i_part_neight && i_rank == i_rank_neight) {
            int i_cell_neight = _prev_cell_cell_extended[3*idx_neight+2];

            int i_unique = _unique_prev_cell_cell_extended[idx_neight]; // Pas la cellule le tableau est unique !!!!!
            // int is_treat = tag_prev_cell_cell_extended    [i_unique  ];

            // printf("i_cell_neight = %i | tag_prev_cell_cell_extended[%i] = %i \n", i_cell_neight, i_unique, is_treat);

            /* From interior */
            // _cell_cell_extended_n[i_cell] += _prev_cell_cell_extended_n[i_cell_neight];

            for(int idx_neight2 = _prev_cell_cell_extended_idx[i_cell_neight]; idx_neight2 < _prev_cell_cell_extended_idx[i_cell_neight+1]; ++idx_neight2) {

              int i_unique2 = _unique_prev_cell_cell_extended[idx_neight2]; // Pas la cellule le tableau est unique !!!!!
              // printf("[%i] [neighbor = %i] _unique_prev_cell_cell_extended[%i] = %i \n", i_cell, i_cell_neight, idx_neight2, i_unique2);
              int is_treat2 = tag_prev_cell_cell_extended    [i_unique2  ];
              if(is_treat2 == 0) {
                _cell_cell_extended_n      [i_cell   ] += 1;
                tag_prev_cell_cell_extended[i_unique2]  = 1;
                icell_to_reset[n_unique_loc++] = i_unique2;
              }
            }

            int idx_border_neight = idx_border_cell[i_cell_neight];
            if(idx_border_neight != -1) {
              // Il faut rajouter les voisins aussi
              // _cell_cell_extended_n[i_cell] += _border_cell_cell_extended_n[idx_border_neight];
              for(int idx_neight2 = _border_cell_cell_extended_idx[idx_border_neight]; idx_neight2 < _border_cell_cell_extended_idx[idx_border_neight+1]; ++idx_neight2) {
                int i_unique2 = _unique_order_border_cell_cell_extended[idx_neight2];
                int is_treat2 = tag_border_cell_cell_extended          [i_unique2 ];
                if(is_treat2 == 0) {
                  _cell_cell_extended_n[i_cell] += 1;
                  tag_border_cell_cell_extended[i_unique2]  = 1;
                  icell_border_to_reset[n_border_unique_loc++] = i_unique2;
                }
              }
            }

            /* Rajout du vrai intérieur */
            // _cell_cell_extended_n[i_cell] += _cell_cell_idx[i_cell_neight+1] - _cell_cell_idx[i_cell_neight];

            for(int idx_neight2 = _cell_cell_idx[i_cell_neight]; idx_neight2 < _cell_cell_idx[i_cell_neight+1]; ++idx_neight2 ) {
              int i_cell_neight2 = _cell_cell[idx_neight2];
              if(tag_interior_cell[i_cell_neight2-1] == 0) {

                _cell_cell_extended_n[i_cell] += 1;

                icell_to_reset_interior[n_interior_unique_loc++] = i_cell_neight2-1;
                tag_interior_cell[i_cell_neight2-1] = 1;
              }
            }

            tag_prev_cell_cell_extended[i_unique] = 1;

          } /* End if same part and same proc */
        } /* End loop neighbor */

        // Reset tag
        // log_trace("Count i_cell : %i -> n_unique_loc = %i | n_interior_unique_loc = %i | n_border_unique_loc = %i | _cell_cell_extended_n = %i\n",
        //                 i_cell, n_unique_loc, n_interior_unique_loc, n_border_unique_loc, _cell_cell_extended_n[i_cell]);
        for(int j = 0; j < n_unique_loc; ++j) {
          tag_prev_cell_cell_extended[icell_to_reset[j]] = 0;
        }
        n_unique_loc = 0;
        for(int j = 0; j < n_interior_unique_loc; ++j) {
          tag_interior_cell[icell_to_reset_interior[j]] = 0;
        }
        n_interior_unique_loc = 0;
        for(int j = 0; j < n_border_unique_loc; ++j) {
          tag_border_cell_cell_extended[icell_border_to_reset[j]] = 0;
        }
        n_border_unique_loc = 0;

      } /* End loop border */

      /* Setup idx and reset */
      int max_neight = 0;
      _cell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell+1] = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];
        max_neight = PDM_MAX(max_neight, _cell_cell_extended_n[i_cell]);
        _cell_cell_extended_n[i_cell] = 0;
      }

      /* Allocate */
      part_ext->cell_cell_extended[i_depth][i_part+shift_part] = (int *) malloc( 4 * _cell_cell_extended_idx[n_cell] * sizeof(int));
      int* _cell_cell_extended = part_ext->cell_cell_extended[i_depth][i_part+shift_part];


      /* Let's go - Second pass to fill */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        int n_unique_loc = 0;
        int n_interior_unique_loc = 0;
        int n_border_unique_loc = 0;
        assert(_cell_cell_extended_n[i_cell] == 0); // All are sorted before

        /* From interior - we add the previous rank */
        for(int idx_neight = _prev_cell_cell_extended_idx[i_cell]; idx_neight < _prev_cell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[4*idx_write  ] = _prev_cell_cell_extended [3*idx_neight  ];
          _cell_cell_extended[4*idx_write+1] = _prev_cell_cell_extended [3*idx_neight+1];
          _cell_cell_extended[4*idx_write+2] = _prev_cell_cell_extended [3*idx_neight+2];
          _cell_cell_extended[4*idx_write+3] = _prev_cell_cell_interface[  idx_neight  ];

          int i_unique2 = _unique_prev_cell_cell_extended[idx_neight]; // Pas la cellule le tableau est unique !!!!!
          tag_prev_cell_cell_extended[i_unique2 ] = 1;
          icell_to_reset[n_unique_loc++] = i_unique2;
        }

        /* From border */
        assert(_border_cell_cell_extended_n[idx_cell] > 0);
        for(int idx_neight = _border_cell_cell_extended_idx[idx_cell]; idx_neight < _border_cell_cell_extended_idx[idx_cell+1]; ++idx_neight) {
          int i_unique2 = _unique_order_border_cell_cell_extended[idx_neight];
          int is_treat2 = tag_border_cell_cell_extended          [i_unique2 ];
          if(is_treat2 == 0) {

            int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
            _cell_cell_extended[4*idx_write  ] = _border_cell_cell_extended [3*idx_neight  ];
            _cell_cell_extended[4*idx_write+1] = _border_cell_cell_extended [3*idx_neight+1];
            _cell_cell_extended[4*idx_write+2] = _border_cell_cell_extended [3*idx_neight+2];
            _cell_cell_extended[4*idx_write+3] = _border_cell_cell_interface[  idx_neight  ];

            tag_border_cell_cell_extended[i_unique2]  = 1;
            icell_border_to_reset[n_border_unique_loc++] = i_unique2;
          }
        }

        /* Now we have to extend the interior */
        for(int idx_neight = _prev_cell_cell_extended_idx[i_cell]; idx_neight < _prev_cell_cell_extended_idx[i_cell+1]; ++idx_neight) {
          int i_rank_neight = _prev_cell_cell_extended[3*idx_neight  ];
          int i_part_neight = _prev_cell_cell_extended[3*idx_neight+1];
          /* We add stencil only if it's local */
          if(i_part+shift_part == i_part_neight && i_rank == i_rank_neight) {

            int i_cell_neight = _prev_cell_cell_extended[3*idx_neight+2];
            int i_unique = _unique_prev_cell_cell_extended[idx_neight]; // Pas la cellule le tableau est unique !!!!!

            /* From interior */
            for(int idx_neight2 = _prev_cell_cell_extended_idx[i_cell_neight]; idx_neight2 < _prev_cell_cell_extended_idx[i_cell_neight+1]; ++idx_neight2) {

              int i_unique2 = _unique_prev_cell_cell_extended[idx_neight2]; // Pas la cellule le tableau est unique !!!!!
              int is_treat2 = tag_prev_cell_cell_extended    [i_unique2  ];
              if(is_treat2 == 0) {
                int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
                _cell_cell_extended[4*idx_write  ] = _prev_cell_cell_extended [3*idx_neight2  ];
                _cell_cell_extended[4*idx_write+1] = _prev_cell_cell_extended [3*idx_neight2+1];
                _cell_cell_extended[4*idx_write+2] = _prev_cell_cell_extended [3*idx_neight2+2];
                _cell_cell_extended[4*idx_write+3] = _prev_cell_cell_interface[  idx_neight2  ];

                // Update tag
                tag_prev_cell_cell_extended[i_unique2]  = 1;
                icell_to_reset[n_unique_loc++] = i_unique2;
              }
            }

            int idx_border_neight = idx_border_cell[i_cell_neight];
            if(idx_border_neight != -1) {
              // Il faut rajouter les voisins aussi
              for(int idx_neight2 = _border_cell_cell_extended_idx[idx_border_neight]; idx_neight2 < _border_cell_cell_extended_idx[idx_border_neight+1]; ++idx_neight2) {

                int i_unique2 = _unique_order_border_cell_cell_extended[idx_neight2];
                int is_treat2 = tag_border_cell_cell_extended          [i_unique2 ];
                if(is_treat2 == 0) {
                  int idx_write = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
                  _cell_cell_extended[4*idx_write  ] = _border_cell_cell_extended [3*idx_neight2  ];
                  _cell_cell_extended[4*idx_write+1] = _border_cell_cell_extended [3*idx_neight2+1];
                  _cell_cell_extended[4*idx_write+2] = _border_cell_cell_extended [3*idx_neight2+2];
                  _cell_cell_extended[4*idx_write+3] = _border_cell_cell_interface[  idx_neight2  ];

                  tag_border_cell_cell_extended[i_unique2]  = 1;
                  icell_border_to_reset[n_border_unique_loc++] = i_unique2;

                }
              }
            }

            /* Rajout du vrai intérieur */
            for(int idx_neight2 = _cell_cell_idx[i_cell_neight]; idx_neight2 < _cell_cell_idx[i_cell_neight+1]; ++idx_neight2 ) {
              int i_cell_neight2 = _cell_cell[idx_neight2];

              if(tag_interior_cell[i_cell_neight2-1] == 0) {

                int idx_write =  _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
                _cell_cell_extended[4*idx_write  ] = i_rank;
                _cell_cell_extended[4*idx_write+1] = i_part;
                _cell_cell_extended[4*idx_write+2] = i_cell_neight2-1;
                _cell_cell_extended[4*idx_write+3] = -40000;

                icell_to_reset_interior[n_interior_unique_loc++] = i_cell_neight2-1;
                tag_interior_cell[i_cell_neight2-1] = 1;
              }
            }
            tag_prev_cell_cell_extended[i_unique] = 1;

          } /* End if same part and same proc */

        } /* End loop neighbor */

        // Reset tag for next cell
        // log_trace("Fill i_cell : %i -> n_unique_loc = %i | n_interior_unique_loc = %i | n_border_unique_loc = %i | _cell_cell_extended_n = %i\n",
        //                 i_cell, n_unique_loc, n_interior_unique_loc, n_border_unique_loc, _cell_cell_extended_n[i_cell]);
        for(int j = 0; j < n_unique_loc; ++j) {
          tag_prev_cell_cell_extended[icell_to_reset[j]] = 0;
        }
        for(int j = 0; j < n_interior_unique_loc; ++j) {
          tag_interior_cell[icell_to_reset_interior[j]] = 0;
        }
        for(int j = 0; j < n_border_unique_loc; ++j) {
          tag_border_cell_cell_extended[icell_border_to_reset[j]] = 0;
        }

      } /* End loop border */
      // free(_tag_cell_is_treated);
      free(tag_border_cell_cell_extended);
      free(tag_dist_neighbor_cell);
      free(tag_prev_cell_cell_extended);
      free(tag_interior_cell);
      free(icell_to_reset);
      free(icell_border_to_reset);
      free(icell_to_reset_interior);

      /* The _cell_cell_extended need to be sorted because many entry is duplicated */
      int* order = malloc( max_neight * sizeof(int));
      int* _ncell_cell_extended_n   = malloc( (n_cell    ) * sizeof(int));
      int* _ncell_cell_extended_idx = malloc( (n_cell + 1) * sizeof(int));
      int* _ncell_cell_extended     = malloc( 4 * _cell_cell_extended_idx[n_cell] * sizeof(int));
      _ncell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

        int beg = _cell_cell_extended_idx[i_cell];
        int n_connect = _cell_cell_extended_idx[i_cell+1] - beg;

        // if(n_connect != _cell_cell_extended_n[i_cell]) {
        //   log_trace("i_cell = %i | n_connect = %i | _cell_cell_extended_n = %i \n", i_cell, n_connect, _cell_cell_extended_n[i_cell]);
        // }

        assert(n_connect == _cell_cell_extended_n[i_cell]);

        PDM_order_lnum_s(&_cell_cell_extended[4*beg], 4, order, n_connect);

        // log_trace("i_cell = %i | n_cell = %i | n_connect = %i \n", i_cell, n_cell, n_connect);
        // PDM_log_trace_array_int(order, )
        // printf(" order[%i] = ", i_cell);
        // for(int i = 0; i < n_connect; ++i) {
        //   printf(" %i", order[i]);
        // }
        // printf("\n");

        _ncell_cell_extended_n  [i_cell  ] = 0;
        _ncell_cell_extended_idx[i_cell+1] = _ncell_cell_extended_idx[i_cell];

        // int idx_unique = -1;
        int last_proc  = -1;
        int last_part  = -1;
        int last_elmt  = -1;
        int last_inte  = -40;
        for(int i = 0; i < n_connect; ++i) {
          int old_order = order[i];
          int curr_proc = _cell_cell_extended[4*(beg+old_order)  ];
          int curr_part = _cell_cell_extended[4*(beg+old_order)+1];
          int curr_cell = _cell_cell_extended[4*(beg+old_order)+2];
          int curr_inte = _cell_cell_extended[4*(beg+old_order)+3];
          int is_same = _is_same_quadruplet(last_proc, last_part, last_elmt, last_inte,
                                            curr_proc, curr_part, curr_cell, curr_inte);

          if(is_same == 0){ // N'est pas le meme
            // idx_unique++;
            last_proc = curr_proc;
            last_part = curr_part;
            last_elmt = curr_cell;
            last_inte = curr_inte;

            int beg_write = 4 * _ncell_cell_extended_idx[i_cell+1];
            // printf(" write in = %i  | beg_write = %i | idx_unique = %i\n", beg_write + idx_unique, beg_write, idx_unique);
            _ncell_cell_extended[beg_write  ] = curr_proc;
            _ncell_cell_extended[beg_write+1] = curr_part;
            _ncell_cell_extended[beg_write+2] = curr_cell;
            _ncell_cell_extended[beg_write+3] = curr_inte;

            /* Increment the new counter */
            _ncell_cell_extended_idx[i_cell+1]++;
            _ncell_cell_extended_n  [i_cell  ]++;
          }
        }
      }

      /* Free old ptr and assign the sort one */
      free(part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part]);
      free(part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part]);
      free(part_ext->cell_cell_extended    [i_depth][i_part+shift_part]);

      _ncell_cell_extended  = realloc(_ncell_cell_extended , 4 * _ncell_cell_extended_idx[n_cell] * sizeof(int));

      part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part] = _ncell_cell_extended_idx;
      part_ext->cell_cell_extended_n  [i_depth][i_part+shift_part] = _ncell_cell_extended_n;

      int *_unique_order_cell_cell_extended = NULL;
      int n_unique = _setup_unique_order_quadruplet(n_cell,
                                                    _ncell_cell_extended_idx,
                                                    _ncell_cell_extended,
                                                    &_unique_order_cell_cell_extended);
      part_ext->unique_order_cell_cell_extended  [i_depth][i_part+shift_part] = _unique_order_cell_cell_extended;
      part_ext->n_unique_order_cell_cell_extended[i_depth][i_part+shift_part] = n_unique;


      quadruplet_to_triplet_and_array(_ncell_cell_extended_idx[n_cell],
                                      _ncell_cell_extended,
                                      &part_ext->cell_cell_interface[i_depth][i_part+shift_part],
                                      &part_ext->cell_cell_extended [i_depth][i_part+shift_part]);
      free(_ncell_cell_extended);
      _ncell_cell_extended = part_ext->cell_cell_extended[i_depth][i_part+shift_part];


      /*
       * Last post-treatment -> Remove entty that come twice : One from interior and another from periodicity
       * Manage the case of one periodic is redondant with and existing join between 2 partitions
       */
      int *_post_cell_cell_extended_idx = malloc((n_cell+1) * sizeof(int));
      int *_ncell_cell_interface = part_ext->cell_cell_interface[i_depth][i_part+shift_part];
      _post_cell_cell_extended_idx[0] = 0;
      for(int i = 0; i < n_cell; ++i) {
        _post_cell_cell_extended_idx[i+1] = _post_cell_cell_extended_idx[i];

        int last_proc_opp    = -10000000;
        int last_part_opp    = -10000000;
        int last_entity_opp  = -10000000;
        int last_intrf_opp   = -10000000;
        for(int idx = _ncell_cell_extended_idx[i]; idx < _ncell_cell_extended_idx[i+1]; ++idx) {

          int i_proc_opp   = _ncell_cell_extended    [3*idx  ];
          int i_part_opp   = _ncell_cell_extended    [3*idx+1];
          int i_entity_opp = _ncell_cell_extended    [3*idx+2];
          int i_intrf_opp  = _ncell_cell_interface[idx];

          if(last_proc_opp   == i_proc_opp   &&
             last_part_opp   == i_part_opp   &&
             last_entity_opp == i_entity_opp)
          {
            if(last_intrf_opp != -40000) {
              // On reecrit
              int idx_write = _post_cell_cell_extended_idx[i+1]++;
              _ncell_cell_extended    [3*idx_write  ] = _ncell_cell_extended    [3*idx  ];
              _ncell_cell_extended    [3*idx_write+1] = _ncell_cell_extended    [3*idx+1];
              _ncell_cell_extended    [3*idx_write+2] = _ncell_cell_extended    [3*idx+2];
              _ncell_cell_interface[  idx_write  ] = _ncell_cell_interface[  idx  ];
            }
          } else {
            last_proc_opp   = i_proc_opp;
            last_part_opp   = i_part_opp;
            last_entity_opp = i_entity_opp;
            last_intrf_opp  = i_intrf_opp;

            // On reecrit
            int idx_write = _post_cell_cell_extended_idx[i+1]++;
            _ncell_cell_extended    [3*idx_write  ] = _ncell_cell_extended    [3*idx  ];
            _ncell_cell_extended    [3*idx_write+1] = _ncell_cell_extended    [3*idx+1];
            _ncell_cell_extended    [3*idx_write+2] = _ncell_cell_extended    [3*idx+2];
            _ncell_cell_interface[  idx_write  ] = _ncell_cell_interface[  idx  ];
          }
        }
      }

      free(_ncell_cell_extended_idx);
      _ncell_cell_extended_idx = _post_cell_cell_extended_idx;
      part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part] = _ncell_cell_extended_idx;

      // Recompute part_ext->cell_cell_extended_idx
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _ncell_cell_extended_n[i_cell] = _ncell_cell_extended_idx[i_cell+1] - _ncell_cell_extended_idx[i_cell];
      }


      if(0 == 1) {
        PDM_log_trace_array_int(_ncell_cell_extended_idx, n_cell+1, "_ncell_cell_extended_idx:: ");
        PDM_log_trace_array_int(_ncell_cell_extended_n  , n_cell  , "_ncell_cell_extended_n:: ");
        PDM_log_trace_array_int(_ncell_cell_extended    , 3 * _ncell_cell_extended_idx[n_cell], "_ncell_cell_extended:: ");
        PDM_log_trace_array_int(_ncell_cell_extended    , 3 * _ncell_cell_extended_idx[n_cell], "_ncell_cell_extended:: ");
        PDM_log_trace_array_int(part_ext->cell_cell_interface[i_depth][i_part+shift_part]    , _ncell_cell_extended_idx[n_cell], "_ncell_cell_interface:: ");

        printf(" --------------------------------------------------------------- \n");
        printf("_ncell_cell_extended :: --------------------- \n");
        for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
          if( _ncell_cell_extended_idx[i_cell+1] > _ncell_cell_extended_idx[i_cell]){
            printf("[%i] i_cell -> %i -->  ", i_part, i_cell);
          }
          for(int idx = _ncell_cell_extended_idx[i_cell]; idx < _ncell_cell_extended_idx[i_cell+1]; ++idx) {
            printf("(%i, %i, intrf = %i) ", _ncell_cell_extended[3*idx+1], _ncell_cell_extended[3*idx+2], part_ext->cell_cell_interface[i_depth][i_part+shift_part][idx]);
          }
          if( _ncell_cell_extended_idx[i_cell+1] > _ncell_cell_extended_idx[i_cell]){
            printf("\n");
          }
        }
        printf("_ncell_cell_extended :: --------------------- END \n");
        printf(" --------------------------------------------------------------- \n");
      }
      // exit(1);

      /* Free */
      free(idx_border_cell);
      free(_border_cell_cell_extended_idx);
      free(order);

      /* Free allocated memory in distant neigbor exhange */
      free(next_cell_cell_extended_n[i_part+shift_part]);
      free(next_cell_cell_extended  [i_part+shift_part]);
      free(next_cell_cell_interface [i_part+shift_part]);
      free(_unique_order_border_cell_cell_extended);

    }
    shift_part += part_ext->n_part[i_domain];
  }

  free(next_cell_cell_extended_n);
  free(next_cell_cell_extended);
  free(next_cell_cell_interface);

}



static
void
_compute_first_extended_cell_graph
(
 PDM_part_extension_t *part_ext
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  int shift_part   = 0;
  int shift_part_g = 0;
  int i_depth_cur = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int n_cell = part_ext->n_cell[i_part+shift_part];
      part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part] = (int *) malloc( (n_cell + 1 ) * sizeof(int));
      part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part] = (int *) malloc( (n_cell     ) * sizeof(int));

      int* _dist_neighbor_cell_n         = part_ext->dist_neighbor_cell_n        [i_part+shift_part];
      int* _dist_neighbor_cell_idx       = part_ext->dist_neighbor_cell_idx      [i_part+shift_part];
      int* _dist_neighbor_cell_desc      = part_ext->dist_neighbor_cell_desc     [i_part+shift_part];
      int* _dist_neighbor_cell_interface = part_ext->dist_neighbor_cell_interface[i_part+shift_part];

      /* Uniquement besoin du graph sur les cellules de bords */
      // int n_cell_border = _dist_neighbor_cell_idx[n_cell];
      int n_cell_border = part_ext->n_cell_border[i_part+shift_part];

      int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part];
      int* _cell_cell_extended_n   = part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part];

      int* _cell_cell_idx = part_ext->cell_cell_idx[i_part+shift_part];
      int* _cell_cell     = part_ext->cell_cell    [i_part+shift_part];

      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell] = 0;
        _cell_cell_extended_n  [i_cell] = 0;
      }
      _cell_cell_extended_idx[n_cell] = 0;

      // printf(" n_cell_border = %i \n", n_cell_border);

      /* First pass to count */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];
        assert(_cell_cell_extended_n[i_cell] == 0); // All are sorted before

        /* From interior */
        _cell_cell_extended_n[i_cell] = _cell_cell_idx[i_cell+1] - _cell_cell_idx[i_cell];

        /* From border */
        assert(_dist_neighbor_cell_n[i_cell] > 0);
        _cell_cell_extended_n[i_cell] += _dist_neighbor_cell_n[i_cell];
      }

      _cell_cell_extended_idx[0] = 0;
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        _cell_cell_extended_idx[i_cell+1] = _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell];
        _cell_cell_extended_n[i_cell] = 0;
      }

      part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part] = (int *) malloc( 4 * _cell_cell_extended_idx[n_cell] * sizeof(int));
      int* _cell_cell_extended = part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part];

      /* Second pass to fill */
      for(int idx_cell = 0; idx_cell < n_cell_border; ++idx_cell) {
        int i_cell = part_ext->border_cell_list[i_part+shift_part][idx_cell];

        /* From interior */
        for(int idx_neight = _cell_cell_idx[i_cell]; idx_neight < _cell_cell_idx[i_cell+1]; ++idx_neight ) {
          int i_cell_neight = _cell_cell[idx_neight];
          int idx_write =  _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[4*idx_write  ] = i_rank;
          _cell_cell_extended[4*idx_write+1] = i_part+shift_part_g;
          _cell_cell_extended[4*idx_write+2] = i_cell_neight-1;
          _cell_cell_extended[4*idx_write+3] = -40000;

        }

        /* From border */
        for(int idx_neight = _dist_neighbor_cell_idx[i_cell]; idx_neight < _dist_neighbor_cell_idx[i_cell+1]; ++idx_neight) {
          int i_rank_neight = _dist_neighbor_cell_desc[3*idx_neight  ];
          int i_part_neight = _dist_neighbor_cell_desc[3*idx_neight+1];
          int i_cell_neight = _dist_neighbor_cell_desc[3*idx_neight+2];
          int idx_write =  _cell_cell_extended_idx[i_cell] + _cell_cell_extended_n[i_cell]++;
          _cell_cell_extended[4*idx_write  ] = i_rank_neight;
          _cell_cell_extended[4*idx_write+1] = i_part_neight;
          _cell_cell_extended[4*idx_write+2] = i_cell_neight;
          _cell_cell_extended[4*idx_write+3] = _dist_neighbor_cell_interface[idx_neight];

        }

        // printf("[%i] _cell_cell_extended_n[%i] = %i\n", i_part, i_cell, _cell_cell_extended_n[i_cell]);
      }


      int* _unique_cell_cell_extended_idx = NULL;
      int* _unique_cell_cell_extended_n   = NULL;
      int* _unique_cell_cell_extended     = NULL;
      _unique_quadruplet(n_cell,
                         _cell_cell_extended_idx,
                         _cell_cell_extended,
                         &_unique_cell_cell_extended_idx,
                         &_unique_cell_cell_extended_n,
                         &_unique_cell_cell_extended);

      free(part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part]);
      free(part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part]);
      free(part_ext->cell_cell_extended    [i_depth_cur][i_part+shift_part]);

      part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part] = _unique_cell_cell_extended_idx;
      part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part] = _unique_cell_cell_extended_n;

      /*
       * Setup a unique among all connecetvity to avoid adding same entity multiple times
       */
      int *_unique_order_cell_cell_extended = NULL;
      int n_unique = _setup_unique_order_quadruplet(n_cell,
                                                    _unique_cell_cell_extended_idx,
                                                    _unique_cell_cell_extended,
                                                    &_unique_order_cell_cell_extended);

      int *_unique_cell_cell_interface = NULL;
      quadruplet_to_triplet_and_array(_unique_cell_cell_extended_idx[n_cell],
                                      _unique_cell_cell_extended,
                                      &_unique_cell_cell_interface,
                                      &part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part]);

      part_ext->unique_order_cell_cell_extended  [i_depth_cur][i_part+shift_part] = _unique_order_cell_cell_extended;
      part_ext->n_unique_order_cell_cell_extended[i_depth_cur][i_part+shift_part] = n_unique;
      part_ext->cell_cell_interface              [i_depth_cur][i_part+shift_part] = _unique_cell_cell_interface;

      // free(_unique_cell_cell_interface);
      free(_unique_cell_cell_extended);

      /*
       * Last post-treatment -> Remove entty that come twice : One from interior and another from periodicity
       * Manage the case of one periodic is redondant with and existing join between 2 partitions
       */
      int *_post_cell_cell_extended_idx = malloc((n_cell+1) * sizeof(int));
      int *_tmp_cell_cell_extended      = part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part];
      _post_cell_cell_extended_idx[0] = 0;
      for(int i = 0; i < n_cell; ++i) {
        _post_cell_cell_extended_idx[i+1] = _post_cell_cell_extended_idx[i];

        int last_proc_opp    = -10000000;
        int last_part_opp    = -10000000;
        int last_entity_opp  = -10000000;
        int last_intrf_opp   = -10000000;
        for(int idx = _unique_cell_cell_extended_idx[i]; idx < _unique_cell_cell_extended_idx[i+1]; ++idx) {

          int i_proc_opp   = _tmp_cell_cell_extended    [3*idx  ];
          int i_part_opp   = _tmp_cell_cell_extended    [3*idx+1];
          int i_entity_opp = _tmp_cell_cell_extended    [3*idx+2];
          int i_intrf_opp  = _unique_cell_cell_interface[idx];

          if(last_proc_opp   == i_proc_opp   &&
             last_part_opp   == i_part_opp   &&
             last_entity_opp == i_entity_opp)
          {
            if(last_intrf_opp != -40000) {
              // On reecrit
              int idx_write = _post_cell_cell_extended_idx[i+1]++;
              _tmp_cell_cell_extended    [3*idx_write  ] = _tmp_cell_cell_extended    [3*idx  ];
              _tmp_cell_cell_extended    [3*idx_write+1] = _tmp_cell_cell_extended    [3*idx+1];
              _tmp_cell_cell_extended    [3*idx_write+2] = _tmp_cell_cell_extended    [3*idx+2];
              _unique_cell_cell_interface[  idx_write  ] = _unique_cell_cell_interface[  idx  ];
            }
          } else {
            last_proc_opp   = i_proc_opp;
            last_part_opp   = i_part_opp;
            last_entity_opp = i_entity_opp;
            last_intrf_opp  = i_intrf_opp;

            // On reecrit
            int idx_write = _post_cell_cell_extended_idx[i+1]++;
            _tmp_cell_cell_extended    [3*idx_write  ] = _tmp_cell_cell_extended    [3*idx  ];
            _tmp_cell_cell_extended    [3*idx_write+1] = _tmp_cell_cell_extended    [3*idx+1];
            _tmp_cell_cell_extended    [3*idx_write+2] = _tmp_cell_cell_extended    [3*idx+2];
            _unique_cell_cell_interface[  idx_write  ] = _unique_cell_cell_interface[  idx  ];
          }
        }
      }

      free(part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part]);
      part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part] = _post_cell_cell_extended_idx;

      // Recompute cell_cell_extended_n
      for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
        part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part][i_cell] = _post_cell_cell_extended_idx[i_cell+1] - _post_cell_cell_extended_idx[i_cell];
      }

      // PDM_log_trace_array_int(_unique_order_cell_cell_extended, _unique_cell_cell_extended_idx[n_cell]  , "_unique_order_cell_cell_extended::");
      if(0 == 1) {
        // PDM_log_trace_array_int(_cell_cell_extended_n  , n_cell  , "t_cell_cell_extended_n::");
        // PDM_log_trace_array_int(_cell_cell_extended_idx, n_cell+1, "t_cell_cell_extended_idx::");
        // PDM_log_trace_array_int(_cell_cell_extended, 3 * _cell_cell_extended_idx[n_cell]  , "t_cell_cell_extended::");
        PDM_log_trace_array_int(part_ext->cell_cell_extended_n  [i_depth_cur][i_part+shift_part], n_cell  , "t_cell_cell_extended_n::");
        PDM_log_trace_array_int(part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part], n_cell+1, "t_cell_cell_extended_idx::");
        PDM_log_trace_array_int(part_ext->cell_cell_extended    [i_depth_cur][i_part+shift_part], 3 * _unique_cell_cell_extended_idx[n_cell]  , "t_cell_cell_extended::");
        PDM_log_trace_array_int(part_ext->cell_cell_interface   [i_depth_cur][i_part+shift_part],     _unique_cell_cell_extended_idx[n_cell]  , "t_cell_cell_interface::");
      }

      if( 0 == 1) {
        printf("_compute_first_extended_cell_graph - part_ext->cell_cell_extended :: --------------------- \n");
        for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
          if( part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part][i_cell+1] > part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part][i_cell]){
            printf("[%i] i_cell -> %i -->  ", i_part, i_cell);
          }
          for(int idx = part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part][i_cell]; idx < part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part][i_cell+1]; ++idx) {
            printf("(%i, %i, intrf = %i) ", part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part][3*idx+1],
                                            part_ext->cell_cell_extended[i_depth_cur][i_part+shift_part][3*idx+2],
                                            part_ext->cell_cell_interface[i_depth_cur][i_part+shift_part][idx]);
          }
          if( part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part][i_cell+1] > part_ext->cell_cell_extended_idx[i_depth_cur][i_part+shift_part][i_cell]){
            printf("\n");
          }
        }
        printf("_compute_first_extended_cell_graph - part_ext->cell_cell_extended :: --------------------- END \n");
      }


    }
    shift_part   += part_ext->n_part              [i_domain];
    shift_part_g += part_ext->n_tot_part_by_domain[i_domain];
  }
}


static
void
_prune_cell_cell_extented
(
  PDM_part_extension_t *part_ext,
  int i_depth
)
{
  // printf("_prune_cell_cell_extented \n");
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int* _cell_cell_extended_idx = part_ext->cell_cell_extended_idx[i_depth][i_part+shift_part];
      int* _cell_cell_extended     = part_ext->cell_cell_extended    [i_depth][i_part+shift_part];
      int* _cell_cell_interface    = part_ext->cell_cell_interface   [i_depth][i_part+shift_part];

      int n_cell      = part_ext->parts[i_domain][i_part].n_cell;
      int s_tot       = _cell_cell_extended_idx[n_cell];

      if( 0 == 1) {
        PDM_log_trace_array_int(_cell_cell_extended_idx, n_cell+1, "(in pruned) _cell_cell_extended_idx::");
        PDM_log_trace_array_int(_cell_cell_extended    , 3 * _cell_cell_extended_idx[n_cell], "(in pruned) _cell_cell_extended::");
        PDM_log_trace_array_int(_cell_cell_interface   ,     _cell_cell_extended_idx[n_cell], "(in pruned) _cell_cell_interface::");
      }

      int *_quad_cell_cell_extended = NULL;
      triplet_to_quadruplet(s_tot,
                            _cell_cell_extended,
                            _cell_cell_interface,
                            &_quad_cell_cell_extended);

      int* order = (int * ) malloc( s_tot * sizeof(int));
      PDM_order_lnum_s(_quad_cell_cell_extended, 4, order, s_tot);

      part_ext->cell_cell_extended_pruned    [i_part+shift_part] = (int * ) malloc(  3 * s_tot * sizeof(int));
      part_ext->cell_cell_extended_pruned_idx[i_part+shift_part] = (int * ) malloc( (n_cell+1) * sizeof(int));
      part_ext->cell_cell_interface_pruned   [i_part+shift_part] = (int * ) malloc(      s_tot * sizeof(int));

      int* _cell_cell_extended_pruned     = part_ext->cell_cell_extended_pruned    [i_part+shift_part];
      int* _cell_cell_extended_pruned_idx = part_ext->cell_cell_extended_pruned_idx[i_part+shift_part];
      int* _cell_cell_interface_pruned    = part_ext->cell_cell_interface_pruned   [i_part+shift_part];

      // PDM_log_trace_array_int(_quad_cell_cell_extended   ,   4 *  s_tot, "_quad_cell_cell_extended::");

      int idx_unique = 0;
      int last_proc  = -1;
      int last_part  = -1;
      int last_elmt  = -1;
      int last_inte  = -40000;
      _cell_cell_extended_pruned_idx[0] = 0;
      _cell_cell_extended_pruned_idx[1] = 0;
      for(int i = 0; i < s_tot; ++i) {
        int old_order = order[i];
        int curr_proc = _quad_cell_cell_extended[4*old_order  ];
        int curr_part = _quad_cell_cell_extended[4*old_order+1];
        int curr_cell = _quad_cell_cell_extended[4*old_order+2];
        int curr_inte = _quad_cell_cell_extended[4*old_order+3];
        int is_same = _is_same_quadruplet(last_proc, last_part, last_elmt, last_inte,
                                          curr_proc, curr_part, curr_cell, curr_inte);


        int is_local = (curr_proc == i_rank) && (curr_part == i_part+shift_part) && (curr_inte == -40000);

        // printf("  last = %i / %i / %i / %i \n", last_proc, last_part, last_elmt, last_inte);
        // printf("  cur  = %i / %i / %i / %i \n", curr_proc, curr_part, curr_cell, curr_inte);
        // printf(" is_same = %i | is_local = %i \n", is_same, is_local);

        // On peut également trie les locaux qui ne serve à rien
        // int is_local = (curr_proc == i_rank) && (curr_part == i_part+shift_part);
        if(is_same == 0 && !is_local){ // N'est pas le meme
          // idx_unique++;
          last_proc = curr_proc;
          last_part = curr_part;
          last_elmt = curr_cell;
          last_inte = curr_inte;

          _cell_cell_extended_pruned [3*idx_unique  ] = curr_proc;
          _cell_cell_extended_pruned [3*idx_unique+1] = curr_part;
          _cell_cell_extended_pruned [3*idx_unique+2] = curr_cell;
          _cell_cell_interface_pruned[idx_unique    ] = curr_inte;
          idx_unique++;

          /* Increment the new counter */
          _cell_cell_extended_pruned_idx[1]++;
        } else {

          last_proc = curr_proc;
          last_part = curr_part;
          last_elmt = curr_cell;
          last_inte = curr_inte;

        }
      }
      free(_quad_cell_cell_extended);

      /* On considère que une cellule est connecté a toutes les autres */
      for(int i = 1; i < n_cell; ++i) {
        _cell_cell_extended_pruned_idx[i+1] = _cell_cell_extended_pruned_idx[i];
      }

      part_ext->cell_cell_extended_pruned [i_part+shift_part] = realloc(part_ext->cell_cell_extended_pruned [i_part+shift_part], 3 * idx_unique * sizeof(int));
      part_ext->cell_cell_interface_pruned[i_part+shift_part] = realloc(part_ext->cell_cell_interface_pruned[i_part+shift_part],     idx_unique * sizeof(int));

      if(0 == 1) {
        _cell_cell_extended_pruned     = part_ext->cell_cell_extended_pruned    [i_part+shift_part];
        _cell_cell_interface_pruned     = part_ext->cell_cell_interface_pruned    [i_part+shift_part];
        PDM_log_trace_array_int(_cell_cell_extended_pruned_idx, n_cell+1, "_cell_cell_extended_pruned_idx:: ");
        PDM_log_trace_array_int(_cell_cell_extended_pruned , 3 * _cell_cell_extended_pruned_idx[1], "_cell_cell_extended_pruned:: ");
        PDM_log_trace_array_int(_cell_cell_interface_pruned,     _cell_cell_extended_pruned_idx[1], "_cell_cell_interface_pruned:: ");
      }
      free(order);
      // exit(1);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  // printf("_prune_cell_cell_extented end \n");
}


static
void
_generate_extended_partition_connectivity
(
 int           n_entity1,
 int           n_entity2,
 int          *entity1_entity1_extended_idx,
 int          *entity1_entity1_extended,
 int          *entity1_entity1_interface,
 PDM_g_num_t  *entity2_ln_to_gn,
 int           n_cur_interface_entity2,
 PDM_g_num_t  *opp_interface_and_gnum_entity2,
 int          *cur_interface_entity2,
 int          *cur_interface_sens,
 int          *entity1_entity2_idx,
 int          *entity1_entity2,
 int          *border_gentity1_entity2_n,
 PDM_g_num_t  *border_gentity1_entity2,
 int          *border_part_and_proc_id,
 int          *border_lentity1_entity2,
 int         **entity2_entity2_extended_idx,
 int         **entity2_entity2_extended,
 int         **entity2_entity2_interface,
 PDM_g_num_t **border_entity2_ln_to_gn
)
{
  PDM_UNUSED(entity1_entity1_extended);
  PDM_UNUSED(entity1_entity1_interface);
  PDM_UNUSED(n_cur_interface_entity2);
  PDM_UNUSED(opp_interface_and_gnum_entity2);
  PDM_UNUSED(cur_interface_entity2);

  PDM_g_num_t* gentity1_entity2 = (PDM_g_num_t *) malloc( entity1_entity2_idx[n_entity1] * sizeof(PDM_g_num_t));
  for(int i = 0; i < entity1_entity2_idx[n_entity1]; ++i) {
    gentity1_entity2[i] = entity2_ln_to_gn[PDM_ABS(entity1_entity2[i])-1];
  }

  int n_neight_tot = entity1_entity1_extended_idx[n_entity1];
  int *_border_gentity1_entity2_idx = (int * ) malloc( (entity1_entity1_extended_idx[n_entity1]+1) * sizeof(int) );

  _border_gentity1_entity2_idx[0] = 0;
  int s_tot = 0;
  for(int i = 0; i < n_neight_tot; ++i) {
    s_tot += border_gentity1_entity2_n[i];
    _border_gentity1_entity2_idx[i+1] = _border_gentity1_entity2_idx[i] + border_gentity1_entity2_n[i];
  }

  if(0 == 1) {
    PDM_log_trace_array_int (_border_gentity1_entity2_idx, n_neight_tot+1, "_border_gcell_face_idx   ::");
    PDM_log_trace_array_int (border_gentity1_entity2_n   , n_neight_tot  , "border_gentity1_entity2_n::");
    PDM_log_trace_array_long(border_gentity1_entity2     , s_tot         , "border_gentity1_entity2  ::");
  }


  /*
   * Prepare and order the current entity ln_to_gn
   *   Cause each entity1 can have connectivity of interior or a new entitity2 (from neightborood)
   */
  PDM_g_num_t* _sorted_entity2_ln_to_gn = (PDM_g_num_t * ) malloc( n_entity2 * sizeof(PDM_g_num_t));
  for(int i_entity2 = 0; i_entity2 < n_entity2; ++i_entity2 ) {
    _sorted_entity2_ln_to_gn[i_entity2] = entity2_ln_to_gn[i_entity2];
  }

  int* order                 = (int *) malloc( n_entity2                      * sizeof(int));
  int* order_entity1_entity2 = (int *) malloc( entity1_entity2_idx[n_entity1] * sizeof(int));
  for(int i = 0; i < n_entity2; ++i) {
    order[i] = i;
  }
  for(int i = 0; i < entity1_entity2_idx[n_entity1]; ++i) {
    order_entity1_entity2[i] = i;
  }

  PDM_sort_long(_sorted_entity2_ln_to_gn, order                , n_entity2                     );
  PDM_sort_long(gentity1_entity2        , order_entity1_entity2, entity1_entity2_idx[n_entity1]);
  // abort(); // Il faut trier le cell_face !!!!! --> Permet de prendre le bon signe aprés !

  if(0 == 1) {
    PDM_log_trace_array_long(gentity1_entity2     , entity1_entity2_idx[n_entity1], "gentity1_entity2     ::");
    PDM_log_trace_array_int (order_entity1_entity2, entity1_entity2_idx[n_entity1], "order_entity1_entity2::");
  }

  // Ce qui compte ici c'est le order de la cellule !!!!
  int* border_order = (int * ) malloc( s_tot * sizeof(int));
  for(int i = 0; i < n_neight_tot; ++i) {
    for(int j = _border_gentity1_entity2_idx[i]; j < _border_gentity1_entity2_idx[i+1]; ++j) {
      border_order[j] = j;
    }
  }

  /*
   * New method by doublet
   */
  PDM_g_num_t *_border_entity2_ln_to_gn_and_interface     = (PDM_g_num_t * ) malloc( 2 * s_tot * sizeof(PDM_g_num_t));
  // PDM_g_num_t *init_border_entity2_ln_to_gn_and_interface = (PDM_g_num_t * ) malloc( 2 * s_tot * sizeof(PDM_g_num_t));
  int         *order_tmp                              = (int         * ) malloc(     s_tot * sizeof(int        ));
  int         *border_gentity1_entity2_interface      = (int         * ) malloc(     s_tot * sizeof(int        ));
  for(int i = 0; i < n_neight_tot; ++i) {
    for(int j = _border_gentity1_entity2_idx[i]; j < _border_gentity1_entity2_idx[i+1]; ++j) {
      _border_entity2_ln_to_gn_and_interface    [2*j  ] = PDM_ABS(border_gentity1_entity2  [j]);
      _border_entity2_ln_to_gn_and_interface    [2*j+1] =         entity1_entity1_interface[i];

      // _border_entity2_ln_to_gn_and_interface    [2*j+1] = PDM_ABS(entity1_entity1_interface[i]);

      border_gentity1_entity2_interface     [j    ] =  entity1_entity1_interface[i];
    }
  }

  int n_unique_test = PDM_order_inplace_unique_long(s_tot, 2, _border_entity2_ln_to_gn_and_interface, order_tmp );
  _border_entity2_ln_to_gn_and_interface = realloc(_border_entity2_ln_to_gn_and_interface, 2 * n_unique_test * sizeof(PDM_g_num_t));
  // _border_entity2_ln_to_gn               = realloc(_border_entity2_ln_to_gn              ,     n_unique_test * sizeof(PDM_g_num_t));

  PDM_g_num_t *_border_entity2_ln_to_gn = malloc( n_unique_test * sizeof(PDM_g_num_t));
  // int         *border_entity2_iterface  = malloc( n_unique_test * sizeof(int        ));

  /*
   * Recopy in seperated array
   */
  for(int i = 0; i < n_unique_test; ++i) {
    _border_entity2_ln_to_gn[i] = _border_entity2_ln_to_gn_and_interface[2*i];
    // border_entity2_iterface [i] = _border_entity2_ln_to_gn_and_interface[2*i+1];
  }
  // free(border_entity2_iterface);
  int n_entity2_unique = n_unique_test;

  if(0 == 1) {
    PDM_log_trace_array_long(_border_entity2_ln_to_gn_and_interface, 2 * n_unique_test, "_border_entity2_ln_to_gn_and_interface ::");
  }

  /* Apply sort */
  PDM_order_array(s_tot, sizeof(int), order_tmp, border_order);

  free(order_tmp );


  int* border_entity1_order  = (int *) malloc( s_tot * sizeof(int));
  int* border_entity2_unique = (int *) malloc( s_tot * sizeof(int));

  // Indices of the first unique in original array in order to find out the new
  int* border_entity2_first_unique = (int *) malloc( n_entity2_unique * sizeof(int));

  PDM_g_num_t last  = 0; // Gn is never 0 because start as 1
  int last_inter    = -4000001;
  int idx_unique = -1;
  int idx_first  = 0;
  for(int i = 0; i < s_tot; ++i) {
    int old_order = border_order[i];
    int pos = PDM_binary_search_gap_int(old_order, _border_gentity1_entity2_idx, n_neight_tot+1);
    border_entity1_order[i] = pos;

    if( (last != PDM_ABS(border_gentity1_entity2[old_order]) ) ||
        (last == PDM_ABS(border_gentity1_entity2[old_order]) && last_inter != border_gentity1_entity2_interface[old_order])) {
      idx_unique = i;
      last       = PDM_ABS(border_gentity1_entity2[old_order]);
      // last_inter = PDM_ABS(border_gentity1_entity2_interface[old_order]);
      last_inter = border_gentity1_entity2_interface[old_order];
      border_entity2_first_unique[idx_first++] = i;
    }
    border_entity2_unique[old_order] = idx_unique;

    // printf(" Search idx -> border_order[%i] = %i  --> idx_unique = %i | last = %i \n", i, old_order, idx_unique, (int)last);
    // border_entity1_order[i] = -1;
    // printf(" Associated cell = %i \n", pos);
  }

  if(0 == 1) {
    // printf(" --------------------------------------------  \n");
    // for(int i = 0; i < s_tot; ++i) {
    //   printf("border_entity2_unique[%i] = %i -> %i \n", i, border_entity2_unique[i], border_entity2_unique[border_order[i]]);
    // }
    // printf(" -------------------------------------------- \n");
    PDM_log_trace_array_long(_border_entity2_ln_to_gn    , n_entity2_unique, "_border_entity2_ln_to_gn::");
    PDM_log_trace_array_int (border_order                , s_tot           , "border_order::");
    PDM_log_trace_array_int (border_entity1_order        , s_tot           , "border_entity1_order::");
    PDM_log_trace_array_int (border_entity2_unique       , s_tot           , "border_entity2_unique::");
    PDM_log_trace_array_int (border_entity2_first_unique, n_entity2_unique , "border_entity2_first_unique::");
  }
  // exit(1);

  /* Pour chaque elements on chercher si il est dans les entity2_ln_to_gn
   *   Si ce n'est pas le cas, c'est un nouvelle element, on parcours la liste unique des bords
   *   donc les eléments des bordes sont également unique
   *   En même temps on réalise la construction du nouveau graph d'échanges pour entity2
   */
  *entity2_entity2_extended_idx = (int * ) malloc( (    n_entity2 + 1   ) * sizeof(int));
  *entity2_entity2_extended     = (int * ) malloc( (3 * n_entity2_unique) * sizeof(int));
  *entity2_entity2_interface    = (int * ) malloc( (    n_entity2_unique) * sizeof(int));
  int* _entity2_entity2_extended_idx = *entity2_entity2_extended_idx;
  int* _entity2_entity2_extended     = *entity2_entity2_extended;
  int* _entity2_entity2_interface    = *entity2_entity2_interface;

  PDM_g_num_t *entity2_extended_gnum = (PDM_g_num_t * ) malloc( 2 * n_entity2_unique * sizeof(PDM_g_num_t));
  int n_entity2_extended = 0;
  _entity2_entity2_extended_idx[0] = 0;
  _entity2_entity2_extended_idx[1] = 0;
  int idx_write = 0;

  for(int i_entity2 = 0; i_entity2 < n_entity2_unique; ++i_entity2) {
    // PDM_g_num_t g_entity2 = _border_entity2_ln_to_gn[i_entity2];
    // int pos = PDM_binary_search_long(g_entity2, _sorted_entity2_ln_to_gn, n_entity2);

    PDM_g_num_t g_entity2 = _border_entity2_ln_to_gn_and_interface[2*i_entity2  ];
    int         i_interf  = _border_entity2_ln_to_gn_and_interface[2*i_entity2+1];

    // log_trace("g_entity2 = "PDM_FMT_G_NUM" | i_interf = %i \n", g_entity2, i_interf);

    /*
     * Management du cas ou on a deja le entity2 mais a travers une interface (donc pas le même gnum mais c'est le gnum_opposé)
     */
    if(i_interf != -40000 && n_cur_interface_entity2 > 0) {
      // int         i_interf_abs  = PDM_ABS(i_interf)-1;
      int         i_interf_abs  = i_interf;

      // Recherche dans le tableau d'interface
      PDM_g_num_t search_elmt[2] = {g_entity2, i_interf_abs};

      int pos_interface = PDM_order_binary_search_long(search_elmt, opp_interface_and_gnum_entity2, 2, n_cur_interface_entity2);
      // printf("    i_interf_abs = %i | g_entity2 = %i --> pos_interface = %i \n", i_interf_abs, (int) g_entity2, pos_interface);
      // printf("    i_interf = %i | g_entity2 = %i --> pos_interface = %i \n", i_interf, (int) g_entity2, pos_interface);
      if(pos_interface != -1) {
        continue;
      }
    }

    int pos = PDM_binary_search_long(g_entity2, _sorted_entity2_ln_to_gn, n_entity2);

    // printf(" \t  -----> Search found [g_entity2=%i] in  _sorted_entity2_ln_to_gn --> pos = %i\n", (int)g_entity2, pos);
    if(pos == -1 || i_interf != -40000) {
    // if(pos == -1 || i_interf != 0) {
      entity2_extended_gnum[2*n_entity2_extended  ] = g_entity2;
      entity2_extended_gnum[2*n_entity2_extended+1] = i_interf;
      n_entity2_extended++;
      // printf("\t\t found [%i] = %i\n", i_entity2, pos);

      _entity2_entity2_extended_idx[1]++;

      int pos_first_unique = border_entity2_first_unique[i_entity2];

      int old_order        = border_order           [pos_first_unique];

      PDM_g_num_t gopp_entity2 = border_gentity1_entity2[  old_order];
      int         opp_entity2  = border_lentity1_entity2[  old_order];
      int         opp_proc     = border_part_and_proc_id[2*old_order  ];
      int         opp_part     = border_part_and_proc_id[2*old_order+1];

      // printf(" border_gentity1_entity2[%i] = %i for g_entity2 = %i \n", old_order, border_gentity1_entity2[old_order], g_entity2);

      // printf(" old_entity1_order = %i | old_entity2_order = %i | new_entity2 = %i \n", old_entity1_order, old_entity2_order, new_entity2);
      // printf(" [pos = %i] [opp_proc = %i | opp_part = %i | opp_entity1 = %i | opp_entity2 = %i] | g_entity2 = %i | gopp_entity2 = %i\n", pos, opp_proc, opp_part, opp_entity1, PDM_ABS(opp_entity2), (int)g_entity2, gopp_entity2);
      // printf(" pos_first_unique = %i | old_order = %i | g_num = %i | g_num_check = %i \n",
      //        pos_first_unique, old_order, g_entity2, gopp_entity2);
      assert(g_entity2 == PDM_ABS(gopp_entity2));

      // printf("_entity2_entity2_extended[%i] = [%i/%i/%i] \n", idx_write, opp_proc, opp_part, PDM_ABS(opp_entity2)-1);

      _entity2_entity2_extended [3*idx_write  ] = opp_proc;
      _entity2_entity2_extended [3*idx_write+1] = opp_part;
      _entity2_entity2_extended [3*idx_write+2] = PDM_ABS(opp_entity2)-1;
      // _entity2_entity2_interface[  idx_write  ] = i_interf;
      assert(i_interf == border_gentity1_entity2_interface[old_order]);
      _entity2_entity2_interface[  idx_write  ] = border_gentity1_entity2_interface[old_order];
      idx_write++;
    }

  }
  // printf(" ----- \n");
  free(_border_entity2_ln_to_gn_and_interface);

  for(int i = 1; i < n_entity2; ++i) {
    _entity2_entity2_extended_idx[i+1] = _entity2_entity2_extended_idx[i];
  }
  *entity2_entity2_extended = realloc(*entity2_entity2_extended,  (3 * n_entity2_extended) * sizeof(int));

  if(0 == 1) {
    _entity2_entity2_extended = *entity2_entity2_extended;
    PDM_log_trace_array_long(entity2_extended_gnum       , 2 * n_entity2_extended    , "entity2_extended_gnum::");
    PDM_log_trace_array_int(_entity2_entity2_extended_idx, n_entity2+1           , "_entity2_entity2_extended_idx::");
    PDM_log_trace_array_int(_entity2_entity2_extended    , 3 * n_entity2_extended, "_entity2_entity2_extended::");
  }

  // free(old_to_new_order_border);
  // exit(1);

  /*
   * Reconstruction de la connectivité de bord
   *   On écrase border_lentity1_entity2
   */
  if(0 == 1 && n_cur_interface_entity2 > 0) {
    PDM_log_trace_array_long(opp_interface_and_gnum_entity2, 2 * n_cur_interface_entity2, "opp_interface_and_gnum_entity2 ::");
  }

  int idx = 0;
  int i_entity2_extented = 0;
  for(int i = 0; i < s_tot; ++i) {
    PDM_g_num_t g_entity2 = PDM_ABS(border_gentity1_entity2[i]);

    if(border_gentity1_entity2_interface[i] != -40000 && n_cur_interface_entity2 > 0) {
      // int         i_interf  = PDM_ABS(border_gentity1_entity2_interface[i])-1;

      // Recherche dans le tableau d'interface
      PDM_g_num_t search_elmt[2] = {g_entity2, border_gentity1_entity2_interface[i]};

      int pos_interface = PDM_order_binary_search_long(search_elmt, opp_interface_and_gnum_entity2, 2, n_cur_interface_entity2);

      // log_trace("Search = %i %i --> pos_interface = %i \n", (int )g_entity2 , border_gentity1_entity2_interface[i], pos_interface);

      if(pos_interface != -1) {
        int sgn    = PDM_SIGN(border_lentity1_entity2[i]); // A aller cherche dans le cell_face de depart
        int sens   = cur_interface_sens[pos_interface];
        border_lentity1_entity2[idx++] = sens * sgn * ( cur_interface_entity2[pos_interface] + 1 ); // Car on shift
        // border_lentity1_entity2[idx++] = - sgn * ( cur_interface_entity2[pos_interface] + 1 ); // Car on shift
        i_entity2_extented++;
        continue;
      }
    }

    /* On cherche d'abord dans le bord - face_extended_gnum is sort by construction */
    // int pos = PDM_binary_search_long(g_entity2, entity2_extended_gnum, n_entity2_extended);

    PDM_g_num_t search_elmt[2] = {g_entity2, border_gentity1_entity2_interface[i]};
    int pos = PDM_order_binary_search_long(search_elmt, entity2_extended_gnum, 2, n_entity2_extended);

    if(pos != -1) {
      int sgn    = PDM_SIGN(border_lentity1_entity2[i]); // A aller cherche dans le cell_face de depart

      // On doit chercher le sign de l'entity2 qu'on garde, car elle impose le signe
      // int i_unique  = border_entity2_unique[i];
      // int old_order = border_order[i_unique];
      // int g_sgn     = PDM_SIGN(border_gentity1_entity2[old_order]);
      // PDM_g_num_t g_num_check = border_gentity1_entity2[old_order];

      // printf(" Border face comming for other proc : pos = %i | idx = %i | old_order = %i | g_sgn = %i | g_num_check = %i | g_entity2 = %i\n", pos, idx, old_order, g_sgn, (int)g_num_check, (int)g_entity2);
      // printf(" Rebuild from exterior [%i] with gnum = "PDM_FMT_G_NUM" and pos : %i - new numbering %i \n ", i, g_entity2, pos, ( pos + n_entity2 + 1 ));

      border_lentity1_entity2[idx++] = sgn * ( pos + n_entity2 + 1 ); // Car on shift
      i_entity2_extented++;
    } else {

      int pos_interior2 = PDM_binary_search_long(g_entity2, _sorted_entity2_ln_to_gn, n_entity2);
      // int pos_interior = PDM_binary_search_long(g_entity2, gentity1_entity2, entity1_entity2_idx[n_entity1]);
      // printf(" Border face comming from interior %i - %i \n", pos_interior, idx);
      // printf(" g_entity2 = %i | pos_interior = %i - pos_interior2 = %i \n", (int)g_entity2, pos_interior, pos_interior2);
      // assert(pos_interior  != -1);
      assert(pos_interior2 != -1);

      int sgn    = PDM_SIGN(border_lentity1_entity2[i]); // A aller cherche dans le cell_face de depart
      int old_pos = order[pos_interior2];
      border_lentity1_entity2[idx++] = sgn * ( old_pos + 1 );

      // PDM_g_num_t old_g_num = _sorted_entity2_ln_to_gn[pos_interior2];
      // printf("Cas 1 : border_lentity1_entity2 [%i] | Cas 2  : %i \n", idx-1, sgn * ( old_pos + 1 ));
      // printf("Cas 1 : border_lentity1_entity2 [%i] \n", sgn * old_g_num);

    }
  }
  // exit(1);
  free(border_gentity1_entity2_interface);

  // Mise à jour
  _border_entity2_ln_to_gn = realloc(_border_entity2_ln_to_gn, n_entity2_extended * sizeof(PDM_g_num_t));
  *border_entity2_ln_to_gn = _border_entity2_ln_to_gn;
  for(int i = 0; i < n_entity2_extended; ++i){
    _border_entity2_ln_to_gn[i] = entity2_extended_gnum[2*i];
  }

  /*
   * Free
   */
  free(gentity1_entity2);
  free(order_entity1_entity2);
  free(entity2_extended_gnum);
  free(border_entity1_order);
  free(border_entity2_unique);
  free(border_entity2_first_unique);
  free(border_order);
  free(order);
  free(_sorted_entity2_ln_to_gn);
  // free(_border_entity2_ln_to_gn); //
  free(_border_gentity1_entity2_idx);
}




// static
// void
// _rebuild_connectivity_cell_face_debug
// (
//   PDM_part_extension_t *part_ext
// )
// {
//   // printf("_rebuild_connectivity_cell_face \n");

//   int n_tot_all_domain = 0;
//   int n_part_loc_all_domain = 0;
//   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
//     n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
//     n_part_loc_all_domain += part_ext->n_part[i_domain];
//   }

//   PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
//                                                            n_part_loc_all_domain,
//                                                            part_ext->n_cell,
//                                                            part_ext->cell_cell_extended_pruned_idx,
//                                                            part_ext->cell_cell_extended_pruned    );

//   /* On doit échanger toutes les connectivités en un seul coup et de manière descendante */
//   /* Donc par exemple cell_face + face_vtx
//    * Ou cell_face + face_edge + edge_vtx
//    * La deduction des autres se fait par transitivité local
//    * Il faut également gerer les conditions limites
//    */

//   /* Prepare */
//   int         **cell_face_n   = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
//   PDM_g_num_t **gcell_face    = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
//   int         **lcell_face    = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
//   PDM_g_num_t **cell_ln_to_gn = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
//   // PDM_g_num_t **cell_flags    = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

//   int shift_part = 0;
//   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
//     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

//       int* cell_face_idx =  part_ext->parts[i_domain][i_part].cell_face_idx;
//       int* cell_face     =  part_ext->parts[i_domain][i_part].cell_face;

//       lcell_face[i_part+shift_part] = cell_face;

//       int n_cell      = part_ext->parts[i_domain][i_part].n_cell;
//       // int n_face      = part_ext->parts[i_domain][i_part].n_face;
//       int s_cell_face = cell_face_idx[n_cell];

//       cell_face_n[i_part+shift_part] = (int         *) malloc( n_cell      * sizeof(int        ));
//       gcell_face [i_part+shift_part] = (PDM_g_num_t *) malloc( s_cell_face * sizeof(PDM_g_num_t));

//       PDM_g_num_t* face_ln_to_gn = part_ext->parts[i_domain][i_part].face_ln_to_gn;

//       cell_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].cell_ln_to_gn;

//       for(int i_cell = 0; i_cell < n_cell; ++i_cell) {
//         cell_face_n[i_part+shift_part][i_cell] = cell_face_idx[i_cell+1] - cell_face_idx[i_cell];
//         for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {
//           int sgn    = PDM_SIGN(cell_face[idx_face]);
//           int i_face = PDM_ABS (cell_face[idx_face])-1;
//           gcell_face[i_part+shift_part][idx_face] = sgn * face_ln_to_gn[i_face];
//           // printf("gcell_face[%i][%i] = %i \n", i_part+shift_part, idx_face, i_part);
//         }
//       }

//     }
//     shift_part += part_ext->n_part[i_domain];
//   }

//   /* Exchange */
//   int         **border_gcell_face_n;
//   PDM_g_num_t **border_gcell_face;
//   PDM_distant_neighbor_exch(dn,
//                             sizeof(PDM_g_num_t),
//                             PDM_STRIDE_VAR,
//                             -1,
//                             cell_face_n,
//                  (void **)  gcell_face,
//                            &border_gcell_face_n,
//                 (void ***) &border_gcell_face);

//   int         **border_lcell_face_n;
//   int         **border_lcell_face;
//   PDM_distant_neighbor_exch(dn,
//                             sizeof(int),
//                             PDM_STRIDE_VAR,
//                             -1,
//                             cell_face_n,
//                  (void **)  lcell_face,
//                            &border_lcell_face_n,
//                 (void ***) &border_lcell_face);

//   /* On fait le cell_ln_to_gn par la même occasion */
//   PDM_distant_neighbor_exch(dn,
//                             sizeof(PDM_g_num_t),
//                             PDM_STRIDE_CST,
//                             1,
//                             NULL,
//                  (void **)  cell_ln_to_gn,
//                             NULL,
//                 (void ***) &part_ext->border_cell_ln_to_gn);

//   free(lcell_face);
//   free(cell_ln_to_gn);

//   /* Post treatment */
//   shift_part = 0;
//   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
//     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

//       int n_cell        = part_ext->n_cell[i_part+shift_part];
//       int n_face        = part_ext->parts[i_domain][i_part].n_face;

//       int* face_face_extended_pruned_idx;
//       int* face_face_extended_pruned;
//       PDM_g_num_t* face_ln_to_gn;
//       _generate_extended_partition_connectivity(n_cell,
//                                                 n_face,
//                                                 part_ext->cell_cell_extended_pruned_idx[i_part+shift_part],
//                                                 part_ext->cell_cell_extended_pruned    [i_part+shift_part],
//                                                 part_ext->parts[i_domain][i_part].face_ln_to_gn,
//                                                 part_ext->parts[i_domain][i_part].cell_face_idx,
//                                                 part_ext->parts[i_domain][i_part].cell_face,
//                                                 border_gcell_face_n    [i_part+shift_part],
//                                                 border_gcell_face      [i_part+shift_part],
//                                                 border_lcell_face      [i_part+shift_part],
//                                                 border_part_and_proc_id[i_part+shift_part],
//                                                 NULL,
//                                                &face_face_extended_pruned_idx,
//                                                &face_face_extended_pruned,
//                                                &face_ln_to_gn);
//       free(face_face_extended_pruned_idx);
//       free(face_face_extended_pruned);
//       free(face_ln_to_gn);
//     }
//     shift_part += part_ext->n_part[i_domain];
//   }

//   /* Pour les faces group on peut faire aussi le gnum location --> Marche pas en multidomain (ou il faut shifter )*/
//   shift_part = 0;
//   for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
//     for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
//       free(cell_face_n[i_part+shift_part]);
//       free(gcell_face[i_part+shift_part]);
//       free(border_gcell_face_n[i_part+shift_part]);
//       free(border_lcell_face_n[i_part+shift_part]);
//       free(border_lcell_face[i_part+shift_part]);
//       free(border_gcell_face[i_part+shift_part]);
//     }
//     shift_part += part_ext->n_part[i_domain];
//   }

//   PDM_distant_neighbor_free(dn);
//   free(cell_face_n);
//   free(gcell_face);
//   free(border_gcell_face_n);
//   free(border_gcell_face);
//   free(border_lcell_face_n);
//   free(border_lcell_face);
//   // printf("_rebuild_connectivity end \n");
// }


static
void
_rebuild_connectivity
(
  PDM_part_extension_t  *part_ext,
  int                   *n_entity1,
  int                   *n_entity2,
  int                  **entity1_entity1_extended_idx,
  int                  **entity1_entity1_extended,
  int                  **entity1_entity1_interface,
  int                  **entity1_entity2_idx,
  int                  **entity1_entity2,
  PDM_g_num_t          **entity2_ln_to_gn,
  int                   *n_cur_interface_entity2,
  PDM_g_num_t          **opp_interface_and_gnum_entity2,
  int                  **cur_interface_entity2,
  int                  **cur_interface_sens,
  int                 ***border_lentity1_entity2_idx,
  int                 ***border_lentity1_entity2,
  int                 ***entity2_entity2_extended_idx,
  int                 ***entity2_entity2_extended,
  int                 ***entity2_entity2_interface,
  PDM_g_num_t         ***border_entity2_ln_to_gn
)
{
  // printf("_rebuild_connectivity \n");

  int i_rank = -1;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);

  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  *entity2_entity2_extended_idx = (int         ** ) malloc( n_part_loc_all_domain * sizeof(int         *));
  *entity2_entity2_extended     = (int         ** ) malloc( n_part_loc_all_domain * sizeof(int         *));
  *entity2_entity2_interface    = (int         ** ) malloc( n_part_loc_all_domain * sizeof(int         *));
  *border_entity2_ln_to_gn      = (PDM_g_num_t ** ) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
  *border_lentity1_entity2_idx  = (int         ** ) malloc( n_part_loc_all_domain * sizeof(int         *));

  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           n_entity1,
                                                           entity1_entity1_extended_idx,
                                                           entity1_entity1_extended);

  /* Prepare */
  int         **entity1_entity2_n = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **gentity1_entity2  = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int         **part_and_proc_id  = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int* _pentity1_entity2_idx = entity1_entity2_idx[i_part+shift_part];
      int* _pentity1_entity2     = entity1_entity2    [i_part+shift_part];

      int pn_entity1        = n_entity1[i_part+shift_part];
      int s_entity1_entity2 = _pentity1_entity2_idx[pn_entity1];

      entity1_entity2_n[i_part+shift_part] = (int         *) malloc( pn_entity1        * sizeof(int        ));
      gentity1_entity2 [i_part+shift_part] = (PDM_g_num_t *) malloc( s_entity1_entity2 * sizeof(PDM_g_num_t));

      part_and_proc_id [i_part+shift_part] = (int         *) malloc( 2 * s_entity1_entity2 * sizeof(int));

      PDM_g_num_t* _pentity2_ln_to_gn = entity2_ln_to_gn[i_part+shift_part];

      for(int i1 = 0; i1 < pn_entity1; ++i1) {
        entity1_entity2_n[i_part+shift_part][i1] = _pentity1_entity2_idx[i1+1] - _pentity1_entity2_idx[i1];
        for(int idx_entity2 = _pentity1_entity2_idx[i1]; idx_entity2 < _pentity1_entity2_idx[i1+1]; ++idx_entity2) {
          int sgn    = PDM_SIGN(_pentity1_entity2[idx_entity2]);
          int i_face = PDM_ABS (_pentity1_entity2[idx_entity2])-1;
          gentity1_entity2[i_part+shift_part][idx_entity2] = sgn * _pentity2_ln_to_gn[i_face];
          part_and_proc_id[i_part+shift_part][2*idx_entity2  ] = i_rank;
          part_and_proc_id[i_part+shift_part][2*idx_entity2+1] = i_part+shift_part;
        }
      }

    }
    shift_part += part_ext->n_part[i_domain];
  }

  /* Exchange */
  int         **border_gentity1_entity2_n;
  PDM_g_num_t **border_gentity1_entity2;
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            entity1_entity2_n,
                 (void **)  gentity1_entity2,
                           &border_gentity1_entity2_n,
                (void ***) &border_gentity1_entity2);

  int         **border_lentity1_entity2_n;
  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            entity1_entity2_n,
                 (void **)  entity1_entity2,
                           &border_lentity1_entity2_n,
                (void ***)  border_lentity1_entity2);

  int         **border_part_and_proc_id;
  int         **border_lentity1_entity2_n_tmp;
  PDM_distant_neighbor_exch(dn,
                            2 * sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            entity1_entity2_n,
                 (void **)  part_and_proc_id,
                           &border_lentity1_entity2_n_tmp,
                (void ***) &border_part_and_proc_id);

  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(border_lentity1_entity2_n_tmp[i_part+shift_part]);
      free(part_and_proc_id[i_part+shift_part]);
    }
    shift_part += part_ext->n_part[i_domain];
  }
  free(border_lentity1_entity2_n_tmp);
  free(part_and_proc_id);

  /* Post treatment */
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      int pn_entity1 = n_entity1[i_part+shift_part];
      int pn_entity2 = n_entity2[i_part+shift_part];

      int         *pborder_lentity1_entity2      = (*border_lentity1_entity2)[i_part+shift_part];
      int         *pentity2_entity2_extended_idx = NULL;
      int         *pentity2_entity2_extended     = NULL;
      int         *pentity2_entity2_interface    = NULL;
      PDM_g_num_t *pentity2_ln_to_gn             = NULL;
      _generate_extended_partition_connectivity(pn_entity1,
                                                pn_entity2,
                                                entity1_entity1_extended_idx  [i_part+shift_part],
                                                entity1_entity1_extended      [i_part+shift_part],
                                                entity1_entity1_interface     [i_part+shift_part],
                                                entity2_ln_to_gn              [i_part+shift_part],
                                                n_cur_interface_entity2       [i_part+shift_part],
                                                opp_interface_and_gnum_entity2[i_part+shift_part],
                                                cur_interface_entity2         [i_part+shift_part],
                                                cur_interface_sens            [i_part+shift_part],
                                                entity1_entity2_idx           [i_part+shift_part],
                                                entity1_entity2               [i_part+shift_part],
                                                border_gentity1_entity2_n     [i_part+shift_part],
                                                border_gentity1_entity2       [i_part+shift_part],
                                                border_part_and_proc_id       [i_part+shift_part],
                                                pborder_lentity1_entity2,
                                               &pentity2_entity2_extended_idx,
                                               &pentity2_entity2_extended,
                                               &pentity2_entity2_interface,
                                               &pentity2_ln_to_gn);

      (*entity2_entity2_extended_idx)[i_part+shift_part] = pentity2_entity2_extended_idx;
      (*entity2_entity2_extended    )[i_part+shift_part] = pentity2_entity2_extended;
      (*entity2_entity2_interface   )[i_part+shift_part] = pentity2_entity2_interface;
      (*border_entity2_ln_to_gn     )[i_part+shift_part] = pentity2_ln_to_gn;

      /* On refait l'index */
      int n_neight_tot = entity1_entity1_extended_idx[i_part+shift_part][pn_entity1];
      (*border_lentity1_entity2_idx)[i_part+shift_part] = PDM_array_new_idx_from_sizes_int(border_lentity1_entity2_n[i_part+shift_part], n_neight_tot);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  /* Pour les faces group on peut faire aussi le gnum location --> Marche pas en multidomain (ou il faut shifter )*/
  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(entity1_entity2_n[i_part+shift_part]);
      free(gentity1_entity2[i_part+shift_part]);
      free(border_gentity1_entity2_n[i_part+shift_part]);
      free(border_lentity1_entity2_n[i_part+shift_part]);
      // free(border_lentity1_entity2[i_part+shift_part]);
      free(border_gentity1_entity2[i_part+shift_part]);
      free(border_part_and_proc_id[i_part+shift_part]);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  PDM_distant_neighbor_free(dn);
  free(entity1_entity2_n);
  free(gentity1_entity2);
  free(border_gentity1_entity2_n);
  free(border_gentity1_entity2);
  free(border_lentity1_entity2_n);
  free(border_part_and_proc_id);
  // printf("_rebuild_connectivity end \n");

}



static
void
_rebuild_connectivity_cell_face
(
  PDM_part_extension_t *part_ext
)
{

  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_cell        = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int          *n_face        = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **cell_face_idx = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **cell_face     = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **face_ln_to_gn = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
  PDM_g_num_t **cell_ln_to_gn = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_cell       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_cell;
      n_face       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
      cell_face_idx[i_part+shift_part] = part_ext->parts[i_domain][i_part].cell_face_idx;
      cell_face    [i_part+shift_part] = part_ext->parts[i_domain][i_part].cell_face;
      face_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
      cell_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].cell_ln_to_gn;
    }
    shift_part += part_ext->n_part[i_domain];
  }


  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           part_ext->n_cell,
                                                           part_ext->cell_cell_extended_pruned_idx,
                                                           part_ext->cell_cell_extended_pruned    );
  assert(part_ext->border_cell_ln_to_gn == NULL);
  /* On fait le cell_ln_to_gn par la même occasion */
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            1,
                            NULL,
                 (void **)  cell_ln_to_gn,
                            NULL,
                (void ***) &part_ext->border_cell_ln_to_gn);

  /*
   * Ici on doit faire intervenir le tableau d'interface
   *   - On doit recevoir plusieurs fois les ln_to_gn qui viennent de différente interface
   */
  if(1 == 0) {
    shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        PDM_log_trace_array_long(part_ext->border_cell_ln_to_gn[shift_part+i_part]      , part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][part_ext->n_cell[shift_part+i_part]], "part_ext->border_cell_ln_to_gn : ");
        PDM_log_trace_array_int (part_ext->cell_cell_interface_pruned[shift_part+i_part], part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][part_ext->n_cell[shift_part+i_part]], "part_ext->cell_cell_interface_pruned : ");
      }
      shift_part += part_ext->n_part[i_domain];
    }
  }


  PDM_distant_neighbor_free(dn);

  // int **border_lcell_face_idx;
  // int **border_lcell_face;
  // int **face_face_extended_idx;
  // int **face_face_extended;
  assert(part_ext->face_face_extended_idx == NULL);
  assert(part_ext->face_face_extended     == NULL);
  _rebuild_connectivity(part_ext,
                        n_cell,
                        n_face,
                        part_ext->cell_cell_extended_pruned_idx,
                        part_ext->cell_cell_extended_pruned,
                        part_ext->cell_cell_interface_pruned,
                        cell_face_idx,
                        cell_face,
                        face_ln_to_gn,
                        part_ext->n_cur_interface_face,
                        part_ext->opp_interface_and_gnum_face,
                        part_ext->cur_interface_face,
                        part_ext->cur_sens_face,
                       &part_ext->border_cell_face_idx,
                       &part_ext->border_cell_face,
                       &part_ext->face_face_extended_idx,
                       &part_ext->face_face_extended,
                       &part_ext->face_face_interface,
                       &part_ext->border_face_ln_to_gn);

  free(n_cell       );
  free(n_face       );
  free(cell_face_idx);
  free(cell_face    );
  free(face_ln_to_gn);
  free(cell_ln_to_gn);
  // exit(1);
}



static
void
_rebuild_connectivity_face_vtx
(
  PDM_part_extension_t *part_ext
)
{

  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_face        = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int          *n_vtx         = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **face_vtx_idx  = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **face_vtx      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **vtx_ln_to_gn  = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_face       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
      n_vtx        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
      face_vtx_idx [i_part+shift_part] = part_ext->parts[i_domain][i_part].face_vtx_idx;
      face_vtx     [i_part+shift_part] = part_ext->parts[i_domain][i_part].face_vtx;
      vtx_ln_to_gn [i_part+shift_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
    }
    shift_part += part_ext->n_part[i_domain];
  }

  assert(part_ext->vtx_vtx_extended_idx == NULL);
  assert(part_ext->vtx_vtx_extended     == NULL);
  _rebuild_connectivity(part_ext,
                        n_face,
                        n_vtx,
                        part_ext->face_face_extended_idx,
                        part_ext->face_face_extended,
                        part_ext->face_face_interface,
                        face_vtx_idx,
                        face_vtx,
                        vtx_ln_to_gn,
                        part_ext->n_cur_interface_vtx,
                        part_ext->opp_interface_and_gnum_vtx,
                        part_ext->cur_interface_vtx,
                        part_ext->cur_sens_vtx,
                       &part_ext->border_face_vtx_idx,
                       &part_ext->border_face_vtx,
                       &part_ext->vtx_vtx_extended_idx,
                       &part_ext->vtx_vtx_extended,
                       &part_ext->vtx_vtx_interface,
                       &part_ext->border_vtx_ln_to_gn);

  free(n_face      );
  free(n_vtx       );
  free(face_vtx_idx);
  free(face_vtx    );
  free(vtx_ln_to_gn);
}


static
void
_rebuild_connectivity_face_edge
(
  PDM_part_extension_t *part_ext
)
{

  // int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    // n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_face         = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int          *n_edge         = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **face_edge_idx  = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **face_edge      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **edge_ln_to_gn  = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_face        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;
      n_edge        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_edge;
      face_edge_idx [i_part+shift_part] = part_ext->parts[i_domain][i_part].face_edge_idx;
      face_edge     [i_part+shift_part] = part_ext->parts[i_domain][i_part].face_edge;
      edge_ln_to_gn [i_part+shift_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
    }
    shift_part += part_ext->n_part[i_domain];
  }

  assert(part_ext->edge_edge_extended_idx == NULL);
  assert(part_ext->edge_edge_extended     == NULL);
  _rebuild_connectivity(part_ext,
                        n_face,
                        n_edge,
                        part_ext->face_face_extended_idx,
                        part_ext->face_face_extended,
                        part_ext->face_face_interface,
                        face_edge_idx,
                        face_edge,
                        edge_ln_to_gn,
                        part_ext->n_cur_interface_edge,
                        part_ext->opp_interface_and_gnum_edge,
                        part_ext->cur_interface_edge,
                        part_ext->cur_sens_edge,
                       &part_ext->border_face_edge_idx,
                       &part_ext->border_face_edge,
                       &part_ext->edge_edge_extended_idx,
                       &part_ext->edge_edge_extended,
                       &part_ext->edge_edge_interface,
                       &part_ext->border_edge_ln_to_gn);

  if(0 == 1) {
    shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

        int _n_edge          = part_ext->parts[i_domain][i_part].n_edge;
        int _n_face          = part_ext->parts[i_domain][i_part].n_face;
        int _n_face_extended = part_ext->face_face_extended_idx[i_part+shift_part][_n_face];
        // int _n_edge_extended = part_ext->edge_edge_extended_idx[i_part+shift_part][_n_edge];
        int         *extended_face_edge_idx = part_ext->border_face_edge_idx[i_part+shift_part];
        int         *extended_face_edge     = part_ext->border_face_edge    [i_part+shift_part];
        PDM_g_num_t *extended_edge_ln_to_gn = part_ext->border_edge_ln_to_gn[i_part+shift_part];

        for(int i_face = 0; i_face < _n_face; ++i_face) {
          printf("i_face = %i - ", i_face);
          for(int idx_edge = face_edge_idx [i_part+shift_part][i_face]; idx_edge < face_edge_idx [i_part+shift_part][i_face+1]; ++idx_edge) {

            int i_edge = PDM_ABS (face_edge[i_part+shift_part][idx_edge])-1;
            int sgn    = PDM_SIGN(face_edge[i_part+shift_part][idx_edge]);
            PDM_g_num_t g_edge = edge_ln_to_gn [i_part+shift_part][i_edge];

            printf(PDM_FMT_G_NUM" ", sgn * g_edge);
          }
          printf("---- \n");
        }
        for(int i_face = 0; i_face < _n_face_extended; ++i_face) {
          printf("extended i_face = %i - ", i_face);
          for(int idx_edge = extended_face_edge_idx[i_face]; idx_edge < extended_face_edge_idx[i_face+1]; ++idx_edge) {

            int i_edge = PDM_ABS (extended_face_edge[idx_edge])-1;
            int sgn    = PDM_SIGN(extended_face_edge[idx_edge]);
            PDM_g_num_t g_edge = 0;
            if(i_edge < _n_edge) {
              g_edge = edge_ln_to_gn [i_part+shift_part][i_edge];
            } else {
              g_edge = extended_edge_ln_to_gn[i_edge-_n_edge];
            }
            printf(PDM_FMT_G_NUM" ", sgn * g_edge);
          }
          printf("---- \n");
        }
      }
      shift_part += part_ext->n_part[i_domain];
    }
  }

  free(n_face       );
  free(n_edge        );
  free(face_edge_idx);
  free(face_edge    );
  free(edge_ln_to_gn );
}


static
void
_rebuild_connectivity_edge_vtx
(
  PDM_part_extension_t *part_ext
)
{

  // int n_tot_all_domain = 0;
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    // n_tot_all_domain      += part_ext->n_tot_part_by_domain[i_domain];
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_edge        = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int          *n_vtx         = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **edge_vtx_idx  = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **edge_vtx      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **vtx_ln_to_gn  = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_edge      [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_edge;
      n_vtx       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
      edge_vtx_idx[i_part+shift_part] = (int *) malloc((n_edge[i_part+shift_part]+1) * sizeof(int));
      // edge_vtx    [i_part+shift_part] = part_ext->parts[i_domain][i_part].edge_vtx;

      edge_vtx    [i_part+shift_part] = (int *) malloc((2 * n_edge[i_part+shift_part]) * sizeof(int));
      int *_edge_vtx = part_ext->parts[i_domain][i_part].edge_vtx;

      for(int i_edge = 0; i_edge < n_edge      [i_part+shift_part]; ++i_edge) {
        edge_vtx[i_part+shift_part][2*i_edge  ] =  _edge_vtx[2*i_edge  ];
        edge_vtx[i_part+shift_part][2*i_edge+1] = -_edge_vtx[2*i_edge+1];
      }

      vtx_ln_to_gn[i_part+shift_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;

      edge_vtx_idx[i_part+shift_part][0] = 0;
      for(int i_edge = 0; i_edge < n_edge[i_part+shift_part]; ++i_edge ){
        edge_vtx_idx[i_part+shift_part][i_edge+1] = edge_vtx_idx[i_part+shift_part][i_edge] + 2;
      }

    }
    shift_part += part_ext->n_part[i_domain];
  }

  assert(part_ext->vtx_vtx_extended_idx == NULL);
  assert(part_ext->vtx_vtx_extended     == NULL);
  _rebuild_connectivity(part_ext,
                        n_edge,
                        n_vtx,
                        part_ext->edge_edge_extended_idx,
                        part_ext->edge_edge_extended,
                        part_ext->edge_edge_interface,
                        edge_vtx_idx,
                        edge_vtx,
                        vtx_ln_to_gn,
                        part_ext->n_cur_interface_vtx,
                        part_ext->opp_interface_and_gnum_vtx,
                        part_ext->cur_interface_vtx,
                        part_ext->cur_sens_vtx,
                       &part_ext->border_edge_vtx_idx,
                       &part_ext->border_edge_vtx,
                       &part_ext->vtx_vtx_extended_idx,
                       &part_ext->vtx_vtx_extended,
                       &part_ext->vtx_vtx_interface,
                       &part_ext->border_vtx_ln_to_gn);

  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(edge_vtx_idx  [i_part+shift_part]);
      free(edge_vtx      [i_part+shift_part]);

      /* Generic algorithm setup a sign on edge_vtx - We remove it by swap */
      int pn_edge = part_ext->parts[i_domain][i_part].n_edge;
      int pn_edge_extented = part_ext->edge_edge_extended_idx[shift_part+i_part][pn_edge];
      for(int i_edge = 0; i_edge < pn_edge_extented; ++i_edge) {
        int i_vtx1 = part_ext->border_edge_vtx[i_part+shift_part][2*i_edge  ];
        int i_vtx2 = part_ext->border_edge_vtx[i_part+shift_part][2*i_edge+1];
        // printf("i_edge = %i | i_vtx1 = %i | i_vtx2 = %i ",i_edge, i_vtx1, i_vtx2);
        if(i_vtx1 < 0 && i_vtx2 > 0) {
          part_ext->border_edge_vtx[i_part+shift_part][2*i_edge  ] = PDM_ABS(i_vtx2);
          part_ext->border_edge_vtx[i_part+shift_part][2*i_edge+1] = PDM_ABS(i_vtx1);
        } else if (i_vtx2 < 0 && i_vtx1 > 0 ) {
          part_ext->border_edge_vtx[i_part+shift_part][2*i_edge  ] = PDM_ABS(i_vtx1);
          part_ext->border_edge_vtx[i_part+shift_part][2*i_edge+1] = PDM_ABS(i_vtx2);
        } else {
          part_ext->border_edge_vtx[i_part+shift_part][2*i_edge  ] = PDM_ABS(i_vtx1);
          part_ext->border_edge_vtx[i_part+shift_part][2*i_edge+1] = PDM_ABS(i_vtx2);
        }

        // printf(" ----> i_vtx1 = %i | i_vtx2 = %i \n",part_ext->border_edge_vtx[i_part+shift_part][2*i_edge  ], part_ext->border_edge_vtx[i_part+shift_part][2*i_edge+1  ]);
      }

    }
    shift_part += part_ext->n_part[i_domain];
  }
  free(n_edge      );
  free(n_vtx       );
  free(edge_vtx_idx);
  free(edge_vtx    );
  free(vtx_ln_to_gn);
}

static
void
_rebuild_face_group
(
  PDM_part_extension_t *part_ext
)
{
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  /* Cell face */
  int          *n_face              = (int         * ) malloc( n_part_loc_all_domain * sizeof(int          ));
  int         **face_group_idg      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  int         **face_group_n        = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *));
  PDM_g_num_t **face_group_ln_to_gn = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
  PDM_g_num_t **face_ln_to_gn_check = (PDM_g_num_t **) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_face       [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_face;

      int* _pface_group_idx = part_ext->parts[i_domain][i_part].face_bound_idx;
      int* _pface_group     = part_ext->parts[i_domain][i_part].face_bound;

      PDM_g_num_t* _pface_group_ln_to_gn = part_ext->parts[i_domain][i_part].face_group_ln_to_gn;
      PDM_g_num_t* _pface_ln_to_gn       = part_ext->parts[i_domain][i_part].face_ln_to_gn;

      /* Creation d'un champs de face contenant les id de group */
      int n_face_group = part_ext->parts[i_domain][i_part].n_face_group;
      int pn_face      = part_ext->parts[i_domain][i_part].n_face;

      int* face_group_idx = malloc( ( pn_face + 1 ) * sizeof(int));
      face_group_n[i_part+shift_part] = malloc( ( pn_face ) * sizeof(int));

      for(int i_face = 0; i_face < pn_face; ++i_face) {
        face_group_n  [i_part+shift_part][i_face] = 0;
      }

      // printf("n_face_group = %i \n", n_face_group);
      for(int i_group = 0; i_group < n_face_group; ++i_group) {
        // printf("_pface_group_idx[%i] = %i --> %i \n", i_group, _pface_group_idx[i_group], _pface_group_idx[i_group+1]);
        for(int idx_face = _pface_group_idx[i_group]; idx_face < _pface_group_idx[i_group+1]; ++idx_face) {
          int i_face = _pface_group[idx_face];
          // printf("[%i] - iface = %i \n ", i_group, i_face);
          face_group_n[i_part+shift_part][i_face-1]++;
        }
      }

      face_group_idx[0] = 0;
      for(int i_face = 0; i_face < pn_face; ++i_face) {
        // printf(" face_group_n[%i] = %i \n", i_face, face_group_n[i_part+shift_part][i_face]);
        face_group_idx[i_face+1] = face_group_idx[i_face] + face_group_n[i_part+shift_part][i_face];
        face_group_n[i_part+shift_part][i_face] = 0;
      }

      face_group_idg     [i_part+shift_part] = malloc( face_group_idx[pn_face] * sizeof(int        ));
      face_group_ln_to_gn[i_part+shift_part] = malloc( face_group_idx[pn_face] * sizeof(PDM_g_num_t));
      face_ln_to_gn_check[i_part+shift_part] = malloc( face_group_idx[pn_face] * sizeof(PDM_g_num_t));

      int         *_face_group_idg      = face_group_idg     [i_part+shift_part];
      PDM_g_num_t *_face_group_ln_to_gn = face_group_ln_to_gn[i_part+shift_part];
      PDM_g_num_t *_face_ln_to_gn_check = face_ln_to_gn_check[i_part+shift_part];

      // PDM_log_trace_array_long(_pface_ln_to_gn, pn_face, "_pface_ln_to_gn::");

      for(int i_group = 0; i_group < n_face_group; ++i_group) {
        for(int idx_face = _pface_group_idx[i_group]; idx_face < _pface_group_idx[i_group+1]; ++idx_face) {
          int i_face    = _pface_group[idx_face];
          int idx_write = face_group_idx[i_face-1] + face_group_n[i_part+shift_part][i_face-1]++;
          _face_group_idg     [idx_write] = i_group;
          _face_group_ln_to_gn[idx_write] = _pface_group_ln_to_gn[idx_face];
          // printf("[%i] _face_ln_to_gn_check[%i] = %i \n", i_part+shift_part, idx_write, (int)_pface_ln_to_gn[i_face-1]);
          _face_ln_to_gn_check[idx_write] = _pface_ln_to_gn[i_face-1];
        }
      }

      // PDM_log_trace_array_long(_face_group_ln_to_gn, face_group_idx[pn_face], "_face_group_ln_to_gn");
      // PDM_log_trace_array_long(_face_ln_to_gn_check, face_group_idx[pn_face], "_face_ln_to_gn_check");

      free(face_group_idx);
    }
    shift_part += part_ext->n_part[i_domain];
  }


  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           n_face,
                                                           part_ext->face_face_extended_idx,
                                                           part_ext->face_face_extended);

  int** border_face_group_idg_n;
  int** border_face_group_idg;
  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            face_group_n,
                  (void **) face_group_idg,
                           &border_face_group_idg_n,
                 (void ***)&border_face_group_idg);

  int         **border_face_group_ln_to_gn_n;
  PDM_g_num_t **border_face_group_ln_to_gn;
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            face_group_n,
                  (void **) face_group_ln_to_gn,
                           &border_face_group_ln_to_gn_n,
                 (void ***)&border_face_group_ln_to_gn);

  int         **border_face_ln_to_gn_check_n;
  PDM_g_num_t **border_face_ln_to_gn_check;
  PDM_distant_neighbor_exch(dn,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            face_group_n,
                  (void **) face_ln_to_gn_check,
                           &border_face_ln_to_gn_check_n,
                 (void ***)&border_face_ln_to_gn_check);

  PDM_distant_neighbor_free(dn);

  /* Post treatment */
  // TODO MANAGEMENT of multiple domain
  // assert(part_ext->n_domain == 1);

  if(part_ext->n_domain > 1) {
    printf("WARNING : _rebuild_face_group is not managed with n_domain > 1 --> n_domain = %i \n", part_ext->n_domain);
  }

  part_ext->border_face_group_idx      = malloc( n_part_loc_all_domain * sizeof(int         *));
  part_ext->border_face_group          = malloc( n_part_loc_all_domain * sizeof(int         *));
  part_ext->border_face_group_ln_to_gn = malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));

  shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      int pn_face = part_ext->parts[i_domain][i_part].n_face;
      int n_face_border = part_ext->face_face_extended_idx[shift_part+i_part][pn_face];

      // En tout rigeur le nombre de group peut être different entre 2 domaines
      int n_face_group = part_ext->parts[i_domain][i_part].n_face_group;
      part_ext->border_face_group_idx[shift_part+i_part] = malloc( (n_face_group+1) * sizeof(int));
      int *_pborder_face_group_idx = part_ext->border_face_group_idx[shift_part+i_part];

      PDM_array_reset_int(_pborder_face_group_idx, n_face_group+1, 0);

      int idx = 0;
      for(int i = 0; i < n_face_border; ++i) {
        for(int j = 0; j < border_face_group_idg_n[shift_part+i_part][i]; ++j) {
          int i_group = border_face_group_idg[shift_part+i_part][idx++];
          _pborder_face_group_idx[i_group+1]++;
        }
      }

      for(int i_group = 1; i_group < n_face_group; ++i_group) {
        _pborder_face_group_idx[i_group+1] += _pborder_face_group_idx[i_group];
      }

      // printf(" _pborder_face_group_idx[%i] = %i\n", n_face_group, _pborder_face_group_idx[n_face_group]);

      part_ext->border_face_group         [shift_part+i_part] = malloc( _pborder_face_group_idx[n_face_group] * sizeof(int        ));
      part_ext->border_face_group_ln_to_gn[shift_part+i_part] = malloc( _pborder_face_group_idx[n_face_group] * sizeof(PDM_g_num_t));

      int         *_pborder_face_group          = part_ext->border_face_group         [shift_part+i_part];
      PDM_g_num_t *_pborder_face_group_ln_to_gn = part_ext->border_face_group_ln_to_gn[shift_part+i_part];

      int* pborder_face_group_n = PDM_array_zeros_int(n_face_group);

      // PDM_log_trace_array_long(part_ext->border_face_ln_to_gn[shift_part+i_part], n_face_border, "border_face_ln_to_gn::");

      idx = 0;
      for(int i = 0; i < n_face_border; ++i) {
        for(int j = 0; j < border_face_group_idg_n[shift_part+i_part][i]; ++j) {
          int i_group = border_face_group_idg[shift_part+i_part][idx];
          int idx_write = _pborder_face_group_idx[i_group] + pborder_face_group_n[i_group]++;
          // _pborder_face_group         [idx_write] = pn_face+idx; // NON car si on a plusieurs group
          _pborder_face_group         [idx_write] = pn_face+i+1;
          _pborder_face_group_ln_to_gn[idx_write] = border_face_group_ln_to_gn[shift_part+i_part][idx];

          PDM_g_num_t g_num_face = border_face_ln_to_gn_check[shift_part+i_part][idx];
          int pos = PDM_binary_search_long(g_num_face, part_ext->border_face_ln_to_gn[shift_part+i_part], n_face_border);

          if (pos < 0) {
            log_trace("failed to find "PDM_FMT_G_NUM"\n", g_num_face);
          }

          // log_trace("Find "PDM_FMT_G_NUM" pos = %i\n", g_num_face, pos);
          // log_trace(" check_ln_to_gn[%i] = "PDM_FMT_G_NUM" | "PDM_FMT_G_NUM" \n", idx, g_num_face, part_ext->border_face_ln_to_gn[shift_part+i_part][pos] );
          assert(g_num_face == part_ext->border_face_ln_to_gn[shift_part+i_part][pos]);
          idx++;
        }
      }

      if(0 == 1) {
        PDM_log_trace_array_int(_pborder_face_group, _pborder_face_group_idx[n_face_group], "_pborder_face_group :: ");
      }

      free(pborder_face_group_n);
    }
    shift_part += part_ext->n_part[i_domain];
  }

  for(int i = 0; i < n_part_loc_all_domain; ++i) {
    free(border_face_group_idg_n[i]);
    free(border_face_group_idg[i]);
    free(border_face_group_ln_to_gn_n[i]);
    free(border_face_group_ln_to_gn[i]);
    free(border_face_ln_to_gn_check_n[i]);
    free(border_face_ln_to_gn_check[i]);
    free(face_group_n[i]);
    free(face_group_idg[i]);
    free(face_group_ln_to_gn[i]);
    free(face_ln_to_gn_check[i]);
  }
  free(border_face_group_idg_n);
  free(border_face_group_idg);
  free(border_face_group_ln_to_gn_n);
  free(border_face_group_ln_to_gn);
  free(border_face_ln_to_gn_check_n);
  free(border_face_ln_to_gn_check);


  free(n_face      );
  free(face_group_idg);
  free(face_group_n);
  free(face_group_ln_to_gn);
  free(face_ln_to_gn_check);

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a part extension structure
 *
 * \param [in]   comm           Communicator
 *
 */

PDM_part_extension_t*
PDM_part_extension_create
(
 const int                n_domain,
 const int               *n_part,
       PDM_extend_type_t  extend_type,
       int                depth,
 const PDM_MPI_Comm       comm,
 const PDM_ownership_t    owner
)
{
  PDM_part_extension_t *part_ext = (PDM_part_extension_t *) malloc(sizeof(PDM_part_extension_t));

  part_ext->n_domain    = n_domain;
  part_ext->n_part      = malloc( n_domain * sizeof(int)); // Make a copy to avoid pb in cython
  for(int i = 0; i < part_ext->n_domain; ++i) {
    part_ext->n_part[i] = n_part[i];
  }
  part_ext->comm        = comm;
  part_ext->owner       = owner;
  part_ext->extend_type = extend_type;
  part_ext->depth       = depth;

  part_ext->n_part_idx    = (int * ) malloc( (n_domain + 1) * sizeof(int));
  part_ext->n_part_g_idx  = (int * ) malloc( (n_domain + 1) * sizeof(int));
  part_ext->parts = malloc(n_domain * sizeof(_part_t *));


  part_ext->n_part_idx  [0] = 0;
  part_ext->n_part_g_idx[0] = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    part_ext->parts[i_domain] = malloc( n_part[i_domain] * sizeof(_part_t));

    for (int i_part = 0; i_part < n_part[i_domain]; i_part++) {
      part_ext->parts[i_domain][i_domain].n_face_group = 0;
      part_ext->parts[i_domain][i_domain].n_edge_group = 0;
    }

    part_ext->n_part_idx[i_domain+1] = part_ext->n_part_idx[i_domain] + part_ext->n_part[i_domain];

    int n_part_l = n_part[i_domain];
    int n_part_g = -100;
    PDM_MPI_Allreduce(&n_part_l, &n_part_g, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    part_ext->n_part_g_idx[i_domain+1] = part_ext->n_part_idx[i_domain] + n_part_g;

  }

  part_ext->neighbor_idx       = NULL;
  part_ext->neighbor_desc      = NULL;
  part_ext->neighbor_interface = NULL;
  part_ext->n_entity_bound     = NULL;
  part_ext->border_cell_list   = NULL;

  part_ext->dist_neighbor_cell_n              = NULL;
  part_ext->dist_neighbor_cell_idx            = NULL;
  part_ext->dist_neighbor_cell_desc           = NULL;
  part_ext->unique_order_dist_neighbor_cell   = NULL;
  part_ext->n_unique_order_dist_neighbor_cell = NULL;

  part_ext->n_tot_part_by_domain = NULL;

  part_ext->entity_cell_idx     = NULL;
  part_ext->entity_cell         = NULL;
  part_ext->entity_cell_opp_idx = NULL;
  part_ext->entity_cell_opp     = NULL;

  part_ext->opp_interface_and_gnum_face = NULL;
  part_ext->opp_interface_and_gnum_edge = NULL;
  part_ext->opp_interface_and_gnum_vtx  = NULL;
  part_ext->cur_interface_face          = NULL;
  part_ext->cur_interface_edge          = NULL;
  part_ext->cur_interface_vtx           = NULL;
  part_ext->cur_sens_face               = NULL;
  part_ext->cur_sens_edge               = NULL;
  part_ext->cur_sens_vtx                = NULL;
  part_ext->n_cur_interface_face        = NULL;
  part_ext->n_cur_interface_edge        = NULL;
  part_ext->n_cur_interface_vtx         = NULL;

  part_ext->cell_cell_extended_idx            = NULL;
  part_ext->cell_cell_extended_n              = NULL;
  part_ext->cell_cell_extended                = NULL;
  part_ext->cell_cell_interface               = NULL;
  part_ext->unique_order_cell_cell_extended   = NULL;
  part_ext->n_unique_order_cell_cell_extended = NULL;
  part_ext->cell_cell_extended_pruned_idx     = NULL;
  part_ext->cell_cell_extended_pruned         = NULL;
  part_ext->cell_cell_interface_pruned        = NULL;

  part_ext->face_face_extended_idx        = NULL;
  part_ext->face_face_extended            = NULL;
  part_ext->face_face_interface           = NULL;

  part_ext->edge_edge_extended_idx        = NULL;
  part_ext->edge_edge_extended            = NULL;
  part_ext->edge_edge_interface           = NULL;

  part_ext->vtx_vtx_extended_idx          = NULL;
  part_ext->vtx_vtx_extended              = NULL;
  part_ext->vtx_vtx_interface             = NULL;

  part_ext->border_cell_face_idx          = NULL;
  part_ext->border_cell_face              = NULL;

  part_ext->border_face_edge_idx          = NULL;
  part_ext->border_face_edge              = NULL;

  part_ext->border_edge_vtx_idx           = NULL;
  part_ext->border_edge_vtx               = NULL;

  part_ext->border_face_vtx_idx           = NULL;
  part_ext->border_face_vtx               = NULL;

  part_ext->border_face_group_idx         = NULL;
  part_ext->border_face_group             = NULL;

  part_ext->border_vtx                    = NULL;

  part_ext->border_cell_ln_to_gn          = NULL;
  part_ext->border_face_ln_to_gn          = NULL;
  part_ext->border_edge_ln_to_gn          = NULL;
  part_ext->border_vtx_ln_to_gn           = NULL;
  part_ext->border_face_group_ln_to_gn    = NULL;

  part_ext->pdi = NULL;

  part_ext->shift_by_domain_cell       = NULL;
  part_ext->shift_by_domain_face       = NULL;
  part_ext->shift_by_domain_edge       = NULL;
  part_ext->shift_by_domain_vtx        = NULL;
  part_ext->shift_by_domain_face_group = NULL;


  part_ext->n_composed_interface      = 0;
  part_ext->composed_interface_idx    = NULL;
  part_ext->composed_interface        = NULL;
  part_ext->composed_ln_to_gn_sorted  = NULL;

  for (int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; i++) {
    part_ext->has_connectivity[i] = PDM_FALSE;
  }

  return part_ext;
}

/**
 *
 * \brief Set data to perform the partitionned mesh extension
 *
 * \warning Deprecated : use the separate setters instead
 *
 */

void
PDM_part_extension_set_part
(
  PDM_part_extension_t *part_ext,
  int                   i_domain,
  int                   i_part,
  int                   n_cell,
  int                   n_face,
  int                   n_face_part_bound,
  int                   n_face_group,
  int                   n_edge,
  int                   n_vtx,
  int                  *cell_face_idx,
  int                  *cell_face,
  int                  *face_cell,
  int                  *face_edge_idx,
  int                  *face_edge,
  int                  *face_vtx_idx,
  int                  *face_vtx,
  int                  *edge_vtx,
  int                  *face_bound_idx,
  int                  *face_bound,
  int                  *face_join_idx,
  int                  *face_join,
  int                  *face_part_bound_proc_idx,
  int                  *face_part_bound_part_idx,
  int                  *face_part_bound,
  int                  *vtx_part_bound_proc_idx,
  int                  *vtx_part_bound_part_idx,
  int                  *vtx_part_bound,
  PDM_g_num_t          *cell_ln_to_gn,
  PDM_g_num_t          *face_ln_to_gn,
  PDM_g_num_t          *edge_ln_to_gn,
  PDM_g_num_t          *vtx_ln_to_gn,
  PDM_g_num_t          *face_group_ln_to_gn,
  double               *vtx_coord
)
{
  if (edge_vtx != NULL) {
    part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = PDM_TRUE;
  }

  if (face_vtx != NULL) {
    part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX] = PDM_TRUE;
  }

  if (face_edge != NULL) {
    part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
  }

  if (cell_face != NULL) {
    part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
  }

  if (face_cell != NULL) {
    part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_FACE_CELL] = PDM_TRUE;
  }

  part_ext->parts[i_domain][i_part].n_cell            = n_cell;
  part_ext->parts[i_domain][i_part].n_face            = n_face;
  part_ext->parts[i_domain][i_part].n_face_part_bound = n_face_part_bound;
  part_ext->parts[i_domain][i_part].n_face_group      = n_face_group;
  part_ext->parts[i_domain][i_part].n_edge            = n_edge;
  part_ext->parts[i_domain][i_part].n_vtx             = n_vtx;

  part_ext->parts[i_domain][i_part].cell_face_idx = cell_face_idx;
  part_ext->parts[i_domain][i_part].cell_face     = cell_face;
  part_ext->parts[i_domain][i_part].face_cell     = face_cell;

  part_ext->parts[i_domain][i_part].face_edge_idx = face_edge_idx;
  part_ext->parts[i_domain][i_part].face_edge     = face_edge;

  part_ext->parts[i_domain][i_part].face_vtx_idx  = face_vtx_idx;
  part_ext->parts[i_domain][i_part].face_vtx      = face_vtx;

  part_ext->parts[i_domain][i_part].edge_vtx      = edge_vtx;

  part_ext->parts[i_domain][i_part].cell_ln_to_gn = cell_ln_to_gn;
  part_ext->parts[i_domain][i_part].face_ln_to_gn = face_ln_to_gn;
  part_ext->parts[i_domain][i_part].edge_ln_to_gn = edge_ln_to_gn;
  part_ext->parts[i_domain][i_part].vtx_ln_to_gn  = vtx_ln_to_gn;

  part_ext->parts[i_domain][i_part].face_part_bound_proc_idx = face_part_bound_proc_idx;
  part_ext->parts[i_domain][i_part].face_part_bound_part_idx = face_part_bound_part_idx;
  part_ext->parts[i_domain][i_part].face_part_bound          = face_part_bound;

  part_ext->parts[i_domain][i_part].vtx_part_bound_proc_idx = vtx_part_bound_proc_idx;
  part_ext->parts[i_domain][i_part].vtx_part_bound_part_idx = vtx_part_bound_part_idx;
  part_ext->parts[i_domain][i_part].vtx_part_bound          = vtx_part_bound;

  part_ext->parts[i_domain][i_part].face_bound_idx      = face_bound_idx;
  part_ext->parts[i_domain][i_part].face_bound          = face_bound;
  part_ext->parts[i_domain][i_part].face_bound_ln_to_gn = face_group_ln_to_gn;

  part_ext->parts[i_domain][i_part].face_group_ln_to_gn = face_group_ln_to_gn;

  part_ext->parts[i_domain][i_part].face_join_idx       = face_join_idx;
  part_ext->parts[i_domain][i_part].face_join           = face_join;

  part_ext->parts[i_domain][i_part].vtx = vtx_coord;
}



void
PDM_part_extension_part_domain_interface_shared_set
(
  PDM_part_extension_t        *part_ext,
  PDM_part_domain_interface_t *pdi
)
{
  part_ext->pdi = pdi;
}

/**
 *
 * \brief Compute a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_compute
(
  PDM_part_extension_t *part_ext
)
{
  /*
   *  A prevoir : reconstruire les "ghost" sans toute la topologie
   *              par exemple pour l'algébre linéaire
   */
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
  PDM_MPI_Comm_size(part_ext->comm, &n_rank);

  int depth = part_ext->depth;

  part_ext->n_tot_part_by_domain = (int *) malloc( part_ext->n_domain * sizeof(int));
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {

    part_ext->n_tot_part_by_domain[i_domain] = -1;
    int n_part_loc = part_ext->n_part[i_domain];
    PDM_MPI_Allreduce(&n_part_loc, &part_ext->n_tot_part_by_domain[i_domain], 1, PDM_MPI_INT, PDM_MPI_SUM, part_ext->comm);

    n_part_loc_all_domain += part_ext->n_part[i_domain];
  }

  if(1 == 1) {
    _offset_parts_by_domain(part_ext, 1);
  }


  part_ext->entity_cell_n    = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->entity_cell_idx  = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->entity_cell      = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->n_cell           = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ));
  part_ext->n_cell_border    = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ));
  part_ext->border_cell_list = (int **) malloc( n_part_loc_all_domain * sizeof(int *));

  assert(part_ext != NULL);
  // printf(" PDM_part_extension_compute : depth = %i | extend_type = %i \n", depth, part_ext->extend_type);

  // assert(part_ext->extend_type == PDM_EXTEND_FROM_FACE);

  _create_cell_cell_graph(part_ext, part_ext->extend_type);

  part_ext->cell_cell_extended_idx            = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->cell_cell_extended_n              = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->cell_cell_extended                = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->unique_order_cell_cell_extended   = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->cell_cell_interface               = (int *** ) malloc( (depth + 1) * sizeof(int  **));
  part_ext->n_unique_order_cell_cell_extended = (int  ** ) malloc( (depth + 1) * sizeof(int   *));

  for(int i_depth = 0; i_depth < depth; ++i_depth) {
    part_ext->cell_cell_extended_idx           [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->cell_cell_extended_n             [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->cell_cell_extended               [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->unique_order_cell_cell_extended  [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->cell_cell_interface              [i_depth] = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
    part_ext->n_unique_order_cell_cell_extended[i_depth] = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ));
  }

  part_ext->cell_cell_extended_pruned_idx = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->cell_cell_extended_pruned     = (int **) malloc( n_part_loc_all_domain * sizeof(int *));
  part_ext->cell_cell_interface_pruned    = (int **) malloc( n_part_loc_all_domain * sizeof(int *));

  // TODO : vtx_cell
  _create_cell_graph_comm(part_ext);

  /*
   * Warm up domain interface --> Usefull to rebuild connectivity inside domain interface
   */
  _compute_other_domain_interface(part_ext);
  // log_trace("_warm_up_domain_interface(part_ext, PDM_BOUND_TYPE_VTX ) \n");
  _warm_up_domain_interface(part_ext, PDM_BOUND_TYPE_VTX );
  // log_trace("_warm_up_domain_interface(part_ext, PDM_BOUND_TYPE_EDGE ) \n");
  _warm_up_domain_interface(part_ext, PDM_BOUND_TYPE_EDGE);
  // log_trace("_warm_up_domain_interface(part_ext, PDM_BOUND_TYPE_FACE ) \n");
  _warm_up_domain_interface(part_ext, PDM_BOUND_TYPE_FACE);
  // exit(1);

  /*
   * Create for first level the proper graph
   */
  _compute_first_extended_cell_graph(part_ext);


  /*
   *   Step 3 : Compute the graph cell with triplet
   *      -> Now we have all the things ok to exchange direclty on cells
   *      -> In order to build the next depth of ghost cell
   *         we need to prepare an extended graph with triplet
   *      -> With the same distant neigbor we exchange for each depth this graph (containing the information of surrounding cells)
   *   --> Init the distant_neightbor from cell
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(part_ext->comm,
                                                           n_part_loc_all_domain,
                                                           part_ext->n_cell,
                                                           part_ext->dist_neighbor_cell_idx,
                                                           part_ext->dist_neighbor_cell_desc);

  // for(int i_depth = 1; i_depth < depth+1; ++i_depth) {
  for(int i_depth = 1; i_depth < depth; ++i_depth) {
    /* Graph compute = local + distant */
    // printf(" ------------------------------------------------- i_depth = %i \n", i_depth);
    _compute_dual_graph(part_ext, dn, i_depth);
  }

  PDM_distant_neighbor_free(dn);

  /*
   * At this step we have for each level the opposite connected cell
   *   In order to setup all ghost cell we need to deduced all descending connectivity
   *       cell -> face -> edge -> vtx
   * Another step is to compute the ln_to_gn in mulitpart context -> (i_domain, ln_to_gn)
   */

  // _prune_cell_cell_extented(part_ext, depth);
  _prune_cell_cell_extented(part_ext, depth-1);

  // if( 0 == 1) {
  //   _rebuild_connectivity_cell_face_debug(part_ext);
  // }
  if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
    _rebuild_connectivity_cell_face(part_ext);
    if (part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE]) {
      assert(part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
      _rebuild_connectivity_face_edge(part_ext);
      _rebuild_connectivity_edge_vtx (part_ext);
    }
    else {
      assert(part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX]);
      _rebuild_connectivity_face_vtx(part_ext);
    }
  } else if (part_ext->extend_type == PDM_EXTEND_FROM_VTX) {
    _rebuild_connectivity_cell_face(part_ext);
    // exit(1);

    if (part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE]) {
      assert(part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
      _rebuild_connectivity_face_edge(part_ext);
      _rebuild_connectivity_edge_vtx (part_ext);
    }
    else {
      assert(part_ext->has_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX]);
      _rebuild_connectivity_face_vtx(part_ext);
    }
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_compute wrong extend_type \n");
  }

  int     *n_vtx     = (int     * ) malloc( n_part_loc_all_domain * sizeof(int     ));
  int     *part_to_domain  = (int     * ) malloc( n_part_loc_all_domain * sizeof(int     ));
  double **vtx_coord = (double ** ) malloc( n_part_loc_all_domain * sizeof(double *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      n_vtx         [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
      vtx_coord     [i_part+shift_part] = part_ext->parts[i_domain][i_part].vtx;
      part_to_domain[i_part+shift_part] = i_domain;
    }
    shift_part += part_ext->n_part[i_domain];
  }

  if(0 == 1) {
    int **old_to_new_indices = NULL;
    _generate_graph_comm_with_extended(part_ext->comm,
                                       n_part_loc_all_domain,
                                       part_to_domain,
                                       part_ext->pdi,
                                       PDM_BOUND_TYPE_VTX,
                                       part_ext->n_part_idx,
                                       part_ext->n_part_g_idx,
                                       n_vtx,
                                       part_ext->vtx_vtx_extended_idx,
                                       part_ext->vtx_vtx_extended,
                                       part_ext->vtx_vtx_interface,
                                       part_ext->n_interface,
                                       part_ext->n_composed_interface,
                                       part_ext->composed_interface_idx,
                                       part_ext->composed_interface,
                                       part_ext->composed_ln_to_gn_sorted,
                                       part_ext->border_vtx_ln_to_gn,
                                       &old_to_new_indices);

    free(part_to_domain);

    /*
     * Apply order - TODO -> Update connectivity + realloc
     */
    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {

      int n_extented_old = part_ext->vtx_vtx_extended_idx[i_part][n_vtx[i_part]];

      int         *new_vtx_vtx_extended_idx = malloc((n_vtx[i_part]+1)  * sizeof(int        ));
      int         *new_vtx_vtx_extended     = malloc(3 * n_extented_old * sizeof(int        ));
      int         *new_vtx_vtx_interface    = malloc(    n_extented_old * sizeof(int        ));
      PDM_g_num_t *new_border_vtx_ln_to_gn  = malloc(    n_extented_old * sizeof(PDM_g_num_t));

      int n_extented_new = 0;
      for(int i = 0; i < n_extented_old; ++i) {
        n_extented_new = PDM_MAX(n_extented_new, old_to_new_indices[i_part][i]+1);
      }

      for(int i = 0; i < n_extented_old; ++i) {
        int idx_write = old_to_new_indices[i_part][i];
        new_vtx_vtx_extended[3*idx_write  ] = part_ext->vtx_vtx_extended [i_part][3*i  ];
        new_vtx_vtx_extended[3*idx_write+1] = part_ext->vtx_vtx_extended [i_part][3*i+1];
        new_vtx_vtx_extended[3*idx_write+2] = part_ext->vtx_vtx_extended [i_part][3*i+2];
        new_vtx_vtx_interface[idx_write]    = part_ext->vtx_vtx_interface[i_part][i];
        new_border_vtx_ln_to_gn[idx_write]  = part_ext->border_vtx_ln_to_gn[i_part][i];
      }

      new_vtx_vtx_extended_idx[0] = 0;
      for(int i = 0; i < n_vtx[i_part]; ++i) {
        new_vtx_vtx_extended_idx[i+1] = n_extented_new;
      }

      free(part_ext->vtx_vtx_extended_idx[i_part]);
      free(part_ext->vtx_vtx_extended    [i_part]);
      free(part_ext->vtx_vtx_interface   [i_part]);
      free(part_ext->border_vtx_ln_to_gn [i_part]);
      part_ext->vtx_vtx_extended_idx[i_part] = new_vtx_vtx_extended_idx;
      part_ext->vtx_vtx_extended    [i_part] = new_vtx_vtx_extended;
      part_ext->vtx_vtx_interface   [i_part] = new_vtx_vtx_interface;
      part_ext->border_vtx_ln_to_gn [i_part] = new_border_vtx_ln_to_gn;
      // PDM_log_trace_array_int(part_ext->vtx_vtx_extended, 3*)

    }
  }
  free(part_to_domain);


  /* Finalize with vertex */
  PDM_distant_neighbor_t* dn_vtx = PDM_distant_neighbor_create(part_ext->comm,
                                                               n_part_loc_all_domain,
                                                               n_vtx,
                                                               part_ext->vtx_vtx_extended_idx,
                                                               part_ext->vtx_vtx_extended);


  assert(part_ext->border_vtx == NULL);


  PDM_distant_neighbor_exch(dn_vtx,
                            sizeof(double),
                            PDM_STRIDE_CST_INTERLACED,
                            3,
                            NULL,
                 (void **)  vtx_coord,
                            NULL,
                (void ***) &part_ext->border_vtx);

  if(part_ext->pdi_neighbor_idx != NULL) {

    int **pdi_neighbor_n = malloc(n_part_loc_all_domain * sizeof(int *));
    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
      pdi_neighbor_n[i_part] = malloc(n_vtx[i_part] * sizeof(int));
      for(int i = 0; i < n_vtx[i_part]; ++i) {
        pdi_neighbor_n[i_part][i] = part_ext->pdi_neighbor_idx[i_part][i+1] - part_ext->pdi_neighbor_idx[i_part][i];
      }
    }

    int **border_neighor_n = NULL;
    int **border_neighor   = NULL;
    PDM_distant_neighbor_exch(dn_vtx,
                              4 * sizeof(int),
                              PDM_STRIDE_VAR_INTERLACED,
                              -1,
                              pdi_neighbor_n,
                   (void **)  part_ext->pdi_neighbor,
                             &border_neighor_n,
                  (void ***) &border_neighor);


    // int **border_neighor_idx = NULL;
    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {

      int *_border_neighor_n = border_neighor_n[i_part];
      // int *_border_neighor   = border_neighor  [i_part];

      int *_vtx_vtx_extended_idx = part_ext->vtx_vtx_extended_idx[i_part];
      // int *_vtx_vtx_extended     = part_ext->vtx_vtx_extended    [i_part];
      // int *_vtx_vtx_interface    = part_ext->vtx_vtx_interface  [i_part];

      int idx_read = 0;
      for(int i = 0; i < n_vtx[i_part]; ++i) {
        // printf("i_vtx = %i \n", i);
        for(int idx = _vtx_vtx_extended_idx[i]; idx < _vtx_vtx_extended_idx[i+1]; ++idx) {
          // printf("\t  ---> Connected with (%i): %i %i %i %i \n", idx,
          //                                                      _vtx_vtx_extended[3*idx  ],
          //                                                      _vtx_vtx_extended[3*idx+1],
          //                                                      _vtx_vtx_extended[3*idx+2],
          //                                                      _vtx_vtx_interface[idx]);

          for(int k = 0; k < _border_neighor_n[idx]; ++k) {
            // printf("\t\t to merge with = %i %i %i %i \n", _border_neighor[4*idx_read], _border_neighor[4*idx_read+1], _border_neighor[4*idx_read+2], _border_neighor[4*idx_read+3]);
            idx_read++;
          }
        }
      }
    }

    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
      free(pdi_neighbor_n            [i_part]);
      free(part_ext->pdi_neighbor_idx[i_part]);
      free(part_ext->pdi_neighbor    [i_part]);
    }
    free(pdi_neighbor_n            );
    free(part_ext->pdi_neighbor_idx);
    free(part_ext->pdi_neighbor    );


    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
      free(border_neighor_n  [i_part]);
      // free(border_neighor_idx[i_part]);
      free(border_neighor    [i_part]);
    }
    free(border_neighor_n  );
    // free(border_neighor_idx);
    free(border_neighor    );
  }


  // // Debug
  // PDM_g_num_t **gnum_check = (PDM_g_num_t ** ) malloc( n_part_loc_all_domain * sizeof(PDM_g_num_t *));
  // PDM_g_num_t **stride  = (int ** ) malloc( n_part_loc_all_domain * sizeof(int *));

  // shift_part = 0;
  // for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
  //   for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
  //     n_vtx        [i_part+shift_part] = part_ext->parts[i_domain][i_part].n_vtx;
  //     gnum_check   [i_part+shift_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;

  //     stride    [i_part+shift_part] = (int *) malloc( part_ext->parts[i_domain][i_part].n_vtx * sizeof(int ));

  //     for( int i_vtx = 0; i_vtx < part_ext->parts[i_domain][i_part].n_vtx; i_vtx++) {
  //        stride    [i_part+shift_part][i_vtx] = 1;
  //     }


  //     PDM_log_trace_array_long(gnum_check   [i_part+shift_part], part_ext->parts[i_domain][i_part].n_vtx        , "vtx_ln_to_gn::");
  //     PDM_log_trace_array_int(part_ext->vtx_vtx_extended_idx[shift_part+i_part], part_ext->parts[i_domain][i_part].n_vtx        , "vtx_vtx_extended_idx::");
  //     PDM_log_trace_array_int(part_ext->vtx_vtx_extended    [shift_part+i_part], 3 * part_ext->vtx_vtx_extended_idx[shift_part+i_part][part_ext->parts[i_domain][i_part].n_vtx], "vtx_vtx_extended_idx::");
  //   }
  //   shift_part += part_ext->n_part[i_domain];
  // }

  // PDM_g_num_t **pgnum_check;
  // int  **recv_strid;
  // // PDM_distant_neighbor_exch(dn_vtx,
  // //                           sizeof(PDM_g_num_t),
  // //                           PDM_STRIDE_CST,
  // //                           1,
  // //                           NULL,
  // //                (void **)  gnum_check,
  // //                           NULL,
  // //               (void ***) &pgnum_check);
  // PDM_distant_neighbor_exch(dn_vtx,
  //                           sizeof(PDM_g_num_t),
  //                           PDM_STRIDE_VAR,
  //                           -1,
  //                           stride,
  //                (void **)  gnum_check,
  //                           &recv_strid,
  //               (void ***) &pgnum_check);

  // shift_part = 0;
  // for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
  //   for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
  //     for(int i_vtx = 0; i_vtx < part_ext->vtx_vtx_extended_idx[shift_part+i_part][part_ext->parts[i_domain][i_part].n_vtx]; ++i_vtx){
  //       printf(" [%i] pgnum_check[%i] = %i \n", i_part, i_vtx, pgnum_check[i_part][i_vtx]);
  //     }
  //   }
  // }

  int n_interface = 0;
  if(part_ext->pdi != NULL) {
    n_interface = PDM_part_domain_interface_n_interface_get(part_ext->pdi);
  }
  double  **translation_vector = malloc(n_interface * sizeof(double *  ));
  double ***rotation_matrix    = malloc(n_interface * sizeof(double ** ));
  double  **rotation_direction = malloc(n_interface * sizeof(double *  ));
  double  **rotation_center    = malloc(n_interface * sizeof(double *  ));
  double   *rotation_angle     = malloc(n_interface * sizeof(double    ));
  for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
    translation_vector[i_interf] = NULL;
    PDM_part_domain_interface_translation_get(part_ext->pdi, i_interf, &translation_vector[i_interf]);

    rotation_matrix[i_interf] = NULL;
    PDM_part_domain_interface_rotation_get   (part_ext->pdi,
                                              i_interf,
                                              &rotation_direction[i_interf],
                                              &rotation_center   [i_interf],
                                              &rotation_angle    [i_interf]);

    if(rotation_center    [i_interf] != NULL) {
      rotation_matrix[i_interf] = malloc(3 * sizeof(double *));
      for(int k = 0; k < 3; ++k) {
        rotation_matrix[i_interf][k] = malloc(3 * sizeof(double));
      }
    }
  }

  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {

    int n_vtx_extended = part_ext->vtx_vtx_extended_idx[i_part][n_vtx[i_part]];
    double *vtx_coord_extended   = part_ext->border_vtx[i_part];
    int    *border_vtx_interface = part_ext->vtx_vtx_interface[i_part];

    for(int i_vtx = 0; i_vtx < n_vtx_extended; ++i_vtx) {
      if(border_vtx_interface[i_vtx] == -40000) {
        continue;
      }
      int i_interf = PDM_ABS(border_vtx_interface[i_vtx])-1;
      if(i_interf < n_interface){

        // printf("i_interf = %i \n", i_interf);
        if(translation_vector[i_interf] != NULL) {
          for(int k = 0; k < 3; ++k) {
            vtx_coord_extended[3*i_vtx+k] += PDM_SIGN(border_vtx_interface[i_vtx]) * translation_vector[i_interf][k];
          }
        } else if(rotation_direction[i_interf] != NULL) {

          double** rot = rotation_matrix[i_interf];

          double x = vtx_coord_extended[3*i_vtx  ];
          double y = vtx_coord_extended[3*i_vtx+1];
          double z = vtx_coord_extended[3*i_vtx+2];

          double xp = x - rotation_center[i_interf][0];
          double yp = y - rotation_center[i_interf][1];
          double zp = z - rotation_center[i_interf][2];

          // rot[0][0] = 1.;
          // rot[0][1] = 0.;
          // rot[0][2] = 0.;

          // rot[1][0] = 0.;
          // rot[1][1] =  cos(PDM_SIGN(border_vtx_interface[i_vtx]) * rotation_angle[i_interf]);
          // rot[1][2] = -sin(PDM_SIGN(border_vtx_interface[i_vtx]) * rotation_angle[i_interf]);

          // rot[2][0] = 0.;
          // rot[2][1] =  sin(PDM_SIGN(border_vtx_interface[i_vtx]) * rotation_angle[i_interf]);
          // rot[2][2] =  cos(PDM_SIGN(border_vtx_interface[i_vtx]) * rotation_angle[i_interf]);

          double c   = cos(PDM_SIGN(border_vtx_interface[i_vtx]) * rotation_angle[i_interf]);
          double s   = sin(PDM_SIGN(border_vtx_interface[i_vtx]) * rotation_angle[i_interf]);
          double omc = 1.-c;
          double xa  = rotation_direction[i_interf][0];
          double ya  = rotation_direction[i_interf][1];
          double za  = rotation_direction[i_interf][2];

          rotation_matrix[i_interf][0][0] = xa * xa * omc + c;
          rotation_matrix[i_interf][0][1] = xa * ya * omc - za * s;
          rotation_matrix[i_interf][0][2] = xa * za * omc + ya * s;

          rotation_matrix[i_interf][1][0] = xa * ya * omc + za * s;
          rotation_matrix[i_interf][1][1] = ya * ya * omc + c;
          rotation_matrix[i_interf][1][2] = ya * za * omc - xa * s;

          rotation_matrix[i_interf][2][0] = xa * za * omc - ya * s;
          rotation_matrix[i_interf][2][1] = ya * za * omc + xa * s;
          rotation_matrix[i_interf][2][2] = za * za * omc + c;

          for(int k = 0; k < 3; ++k) {
            // vtx_coord_extended[3*i_vtx+k] = PDM_SIGN(border_vtx_interface[i_vtx]) * (rot[k][0]*x + rot[k][1]*y + rot[k][2]*z);
            vtx_coord_extended[3*i_vtx+k] = rotation_center[i_interf][k] + (rot[k][0]*xp + rot[k][1]*yp + rot[k][2]*zp);
          }

        }
      } else {
        // int l_interf = PDM_ABS(border_vtx_interface[i_vtx]) - n_interface - 1;
        int l_interf = PDM_binary_search_long(border_vtx_interface[i_vtx], part_ext->composed_ln_to_gn_sorted, part_ext->n_composed_interface);
        for(int idx_comp = part_ext->composed_interface_idx[l_interf]; idx_comp < part_ext->composed_interface_idx[l_interf+1]; ++idx_comp) {
          int i_tr = part_ext->composed_interface[idx_comp];
          for(int k = 0; k < 3; ++k) {
            vtx_coord_extended[3*i_vtx+k] += PDM_SIGN(i_tr) * translation_vector[PDM_ABS(i_tr)-1][k];
          }
        }
      }
    }
  }


  for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
    if(translation_vector[i_interf] != NULL) {
      free(translation_vector[i_interf]);
    }
    if(rotation_center    [i_interf] != NULL) {
      for(int k = 0; k < 3; ++k) {
        free(rotation_matrix[i_interf][k]);
      }
      free(rotation_matrix[i_interf]);
    }
  }

  free(translation_vector);
  free(rotation_matrix);
  free(rotation_direction);
  free(rotation_center);
  free(rotation_angle);


  free(n_vtx);
  free(vtx_coord);
  // log_trace(" PDM_part_extension_compute flag END \n");

  PDM_distant_neighbor_free(dn_vtx);
  // printf(" PDM_part_extension_compute end \n");

  /* Condition limite - Face uniquement pour l'instant */
  _rebuild_face_group(part_ext);

  /*
   * Debug
   */
  if(0 == 1) {
    shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

        int     pn_vtx        = part_ext->parts[i_domain][i_part].n_vtx;
        int     pn_face       = part_ext->parts[i_domain][i_part].n_face;

        double *pvtx_coord    = part_ext->parts[i_domain][i_part].vtx;
        int    *pface_vtx_idx = part_ext->parts[i_domain][i_part].face_vtx_idx;
        int    *pface_vtx     = part_ext->parts[i_domain][i_part].face_vtx;

        int n_vtx_extended = part_ext->vtx_vtx_extended_idx[shift_part+i_part][pn_vtx];
        int n_face_extended         = part_ext->face_face_extended_idx[shift_part+i_part][pn_face];
        int* pface_vtx_extented_idx = part_ext->border_face_vtx_idx   [shift_part+i_part];
        int* pface_vtx_extented     = part_ext->border_face_vtx       [shift_part+i_part];
        double* pvtx_coord_extented = part_ext->border_vtx            [shift_part+i_part];

        PDM_g_num_t *pvtx_ln_to_gn        = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
        PDM_g_num_t *pface_ln_to_gn       = part_ext->parts[i_domain][i_part].face_ln_to_gn;
        PDM_g_num_t *border_vtx_ln_to_gn  = part_ext->border_vtx_ln_to_gn [i_part+shift_part];
        PDM_g_num_t *border_face_ln_to_gn = part_ext->border_face_ln_to_gn[i_part+shift_part];

        int n_vtx_tot = pn_vtx + n_vtx_extended;
        int n_face_tot = pn_face+n_face_extended;

        PDM_g_num_t *concat_face_ln_to_gn = malloc( (n_face_tot+1) * sizeof(PDM_g_num_t));
        int *concat_face_vtx_idx  = malloc( (n_face_tot+1) * sizeof(int));
        int *concat_face_vtx      = malloc( (pface_vtx_idx[pn_face] + pface_vtx_extented_idx[n_face_extended])* sizeof(int));


        int *color_face = (int *) malloc(n_face_tot * sizeof(int));
        concat_face_vtx_idx[0] = 0;
        for(int i_face = 0; i_face < pn_face; ++i_face) {
          concat_face_vtx_idx[i_face+1] = concat_face_vtx_idx[i_face];
          for(int idx_vtx = pface_vtx_idx[i_face]; idx_vtx < pface_vtx_idx[i_face+1]; ++idx_vtx) {
            concat_face_vtx[concat_face_vtx_idx[i_face+1]++] = pface_vtx[idx_vtx];
          }
          color_face         [i_face] = 0;
          concat_face_ln_to_gn[i_face] = pface_ln_to_gn[i_face];
        }

        for(int i_face = 0; i_face < n_face_extended; ++i_face) {
          int l_face = i_face + pn_face;
          concat_face_vtx_idx[l_face+1] = concat_face_vtx_idx[l_face];
          for(int idx_vtx = pface_vtx_extented_idx[i_face]; idx_vtx < pface_vtx_extented_idx[i_face+1]; ++idx_vtx) {
            concat_face_vtx[concat_face_vtx_idx[l_face+1]++] = PDM_ABS(pface_vtx_extented[idx_vtx]);
          }
          color_face[l_face] = 1;
          concat_face_ln_to_gn[l_face] = border_face_ln_to_gn[i_face];
        }

        PDM_g_num_t *concat_vtx_ln_to_gn  = malloc(    n_vtx_tot  * sizeof(PDM_g_num_t));
        double      *concat_vtx_coord     = malloc(3 * n_vtx_tot  * sizeof(double     ));
        for(int i_vtx = 0; i_vtx < 3 * n_vtx_tot; ++i_vtx) {
          concat_vtx_coord[i_vtx] = -10000.;
        }

        for(int i_vtx = 0; i_vtx < pn_vtx; ++i_vtx) {
          concat_vtx_ln_to_gn[i_vtx] = pvtx_ln_to_gn[i_vtx];
          concat_vtx_coord[3*i_vtx  ] = pvtx_coord[3*i_vtx  ];
          concat_vtx_coord[3*i_vtx+1] = pvtx_coord[3*i_vtx+1];
          concat_vtx_coord[3*i_vtx+2] = pvtx_coord[3*i_vtx+2];
        }

        for(int i_vtx = 0; i_vtx < n_vtx_extended; ++i_vtx) {
          concat_vtx_ln_to_gn[pn_vtx+i_vtx] = border_vtx_ln_to_gn[i_vtx];
          concat_vtx_coord[3*(pn_vtx+i_vtx)  ] = pvtx_coord_extented[3*i_vtx  ];
          concat_vtx_coord[3*(pn_vtx+i_vtx)+1] = pvtx_coord_extented[3*i_vtx+1];
          concat_vtx_coord[3*(pn_vtx+i_vtx)+2] = pvtx_coord_extented[3*i_vtx+2];
        }


        char filename[999];
        sprintf(filename, "out_part_extension_face_vtx_%i_%i_%i.vtk", i_domain, i_part, i_rank);
        PDM_vtk_write_polydata(filename,
                               n_vtx_tot,
                               concat_vtx_coord,
                               concat_vtx_ln_to_gn,
                               n_face_tot,
                               concat_face_vtx_idx,
                               concat_face_vtx,
                               concat_face_ln_to_gn,
                               NULL);
        free(color_face);
        free(concat_face_vtx_idx);
        free(concat_face_vtx);
        free(concat_vtx_coord);
        free(concat_vtx_ln_to_gn);


      }
      shift_part += part_ext->n_part[i_domain];
    }
  }




  if(1 == 1) {
    _offset_parts_by_domain(part_ext, -1);
    _offset_results_by_domain(part_ext);
  }


}


/**
 *
 * \brief Free a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
)
{
  if (part_ext == NULL) {
    return;
  }

  if(part_ext->n_tot_part_by_domain != NULL) {
    free(part_ext->n_tot_part_by_domain);
    part_ext->n_tot_part_by_domain = NULL;
  }

  if(part_ext->neighbor_idx != NULL) {
    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(part_ext->neighbor_idx       [i_part+shift_part]);
        free(part_ext->neighbor_desc      [i_part+shift_part]);
        free(part_ext->neighbor_interface [i_part+shift_part]);

        free(part_ext->dist_neighbor_cell_n           [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_idx         [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_desc        [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_interface   [i_part+shift_part]);
        free(part_ext->unique_order_dist_neighbor_cell[i_part+shift_part]);

        free(part_ext->entity_cell_opp_idx[i_part+shift_part]);
        free(part_ext->entity_cell_opp_n  [i_part+shift_part]);
        free(part_ext->entity_cell_opp    [i_part+shift_part]);

        free(part_ext->border_cell_list    [i_part+shift_part]);

        // if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        free(part_ext->entity_cell_idx[i_part+shift_part]);
        free(part_ext->entity_cell_n  [i_part+shift_part]);
        free(part_ext->entity_cell    [i_part+shift_part]);
        // }


        if(part_ext->opp_interface_and_gnum_face[i_part+shift_part] != NULL) {
          free(part_ext->opp_interface_and_gnum_face[i_part+shift_part]);
        }

        if(part_ext->opp_interface_and_gnum_edge[i_part+shift_part] != NULL) {
          free(part_ext->opp_interface_and_gnum_edge[i_part+shift_part]);
        }

        if(part_ext->opp_interface_and_gnum_vtx[i_part+shift_part] != NULL) {
          free(part_ext->opp_interface_and_gnum_vtx[i_part+shift_part]);
        }

        if(part_ext->cur_interface_face[i_part+shift_part] != NULL) {
          free(part_ext->cur_interface_face[i_part+shift_part]);
        }

        if(part_ext->cur_interface_edge[i_part+shift_part] != NULL) {
          free(part_ext->cur_interface_edge[i_part+shift_part]);
        }

        if(part_ext->cur_interface_vtx[i_part+shift_part] != NULL) {
          free(part_ext->cur_interface_vtx[i_part+shift_part]);
        }

        if(part_ext->cur_sens_face[i_part+shift_part] != NULL) {
          free(part_ext->cur_sens_face[i_part+shift_part]);
        }

        if(part_ext->cur_sens_edge[i_part+shift_part] != NULL) {
          free(part_ext->cur_sens_edge[i_part+shift_part]);
        }

        if(part_ext->cur_sens_vtx[i_part+shift_part] != NULL) {
          free(part_ext->cur_sens_vtx[i_part+shift_part]);
        }


        if(part_ext->face_face_extended_idx != NULL) {
          free(part_ext->face_face_extended_idx[i_part+shift_part]);
        }

        if(part_ext->face_face_extended != NULL) {
          free(part_ext->face_face_extended[i_part+shift_part]);
        }

        if(part_ext->edge_edge_extended_idx != NULL) {
          free(part_ext->edge_edge_extended_idx[i_part+shift_part]);
        }

        if(part_ext->edge_edge_extended != NULL) {
          free(part_ext->edge_edge_extended[i_part+shift_part]);
        }

        if(part_ext->vtx_vtx_extended_idx != NULL) {
          free(part_ext->vtx_vtx_extended_idx[i_part+shift_part]);
        }

        if(part_ext->vtx_vtx_extended != NULL) {
          free(part_ext->vtx_vtx_extended[i_part+shift_part]);
        }

        if(part_ext->owner == PDM_OWNERSHIP_KEEP) {

          if(part_ext->border_cell_face_idx != NULL) {
            free(part_ext->border_cell_face_idx[i_part+shift_part]);
            free(part_ext->border_cell_face    [i_part+shift_part]);
          }

          if(part_ext->border_face_edge_idx != NULL) {
            free(part_ext->border_face_edge_idx [i_part+shift_part]);
            free(part_ext->border_face_edge     [i_part+shift_part]);
          }

          if(part_ext->border_edge_vtx_idx != NULL) {
            free(part_ext->border_edge_vtx_idx[i_part+shift_part]);
            free(part_ext->border_edge_vtx    [i_part+shift_part]);
          }

          if(part_ext->border_face_vtx_idx != NULL) {
            free(part_ext->border_face_vtx_idx[i_part+shift_part]);
            free(part_ext->border_face_vtx    [i_part+shift_part]);
          }

          if(part_ext->border_face_group_idx != NULL) {
            free(part_ext->border_face_group_idx[i_part+shift_part]);
            free(part_ext->border_face_group    [i_part+shift_part]);
          }

          if(part_ext->border_vtx != NULL) {
            free(part_ext->border_vtx[i_part+shift_part]);
          }

          if(part_ext->border_cell_ln_to_gn != NULL) {
            free(part_ext->border_cell_ln_to_gn[i_part+shift_part]);
          }

          if(part_ext->border_face_ln_to_gn != NULL) {
            free(part_ext->border_face_ln_to_gn[i_part+shift_part]);
            free(part_ext->face_face_interface  [i_part+shift_part]);
          }

          if(part_ext->border_edge_ln_to_gn != NULL) {
            free(part_ext->border_edge_ln_to_gn[i_part+shift_part]);
            free(part_ext->edge_edge_interface[i_part+shift_part]);
          }

          if(part_ext->border_vtx_ln_to_gn != NULL) {
            free(part_ext->border_vtx_ln_to_gn[i_part+shift_part]);
            free(part_ext->vtx_vtx_interface  [i_part+shift_part]);
          }

          if(part_ext->border_face_group_ln_to_gn != NULL) {
            free(part_ext->border_face_group_ln_to_gn[i_part+shift_part]);
          }

        }

        free(part_ext->cell_cell_idx[i_part+shift_part]);
        free(part_ext->cell_cell    [i_part+shift_part]);

      }

      shift_part += part_ext->n_part[i_domain];
    }
    free(part_ext->n_entity_bound);
  }
  free(part_ext->neighbor_idx);
  free(part_ext->neighbor_desc);
  free(part_ext->neighbor_interface);
  free(part_ext->n_cell);
  free(part_ext->n_cell_border);
  free(part_ext->border_cell_list);

  free(part_ext->cell_cell_idx);
  free(part_ext->cell_cell);

  if(part_ext->face_face_extended_idx != NULL) {
    free(part_ext->face_face_extended_idx);
    free(part_ext->face_face_extended);
    free(part_ext->face_face_interface);
  }

  if(part_ext->opp_interface_and_gnum_face != NULL) {
    free(part_ext->opp_interface_and_gnum_face);
    free(part_ext->cur_interface_face);
    free(part_ext->cur_sens_face);
    free(part_ext->n_cur_interface_face);
  }
  if(part_ext->opp_interface_and_gnum_edge != NULL) {
    free(part_ext->opp_interface_and_gnum_edge);
    free(part_ext->cur_interface_edge);
    free(part_ext->cur_sens_edge);
    free(part_ext->n_cur_interface_edge);
  }
  if(part_ext->opp_interface_and_gnum_vtx != NULL) {
    free(part_ext->opp_interface_and_gnum_vtx);
    free(part_ext->cur_interface_vtx);
    free(part_ext->cur_sens_vtx);
    free(part_ext->n_cur_interface_vtx);
  }


  if(part_ext->edge_edge_extended_idx != NULL) {
    free(part_ext->edge_edge_extended_idx);
    free(part_ext->edge_edge_extended);
    free(part_ext->edge_edge_interface);
  }

  if(part_ext->vtx_vtx_extended_idx != NULL) {
    free(part_ext->vtx_vtx_extended_idx);
    free(part_ext->vtx_vtx_extended);
    free(part_ext->vtx_vtx_interface);
  }

  /* Peu import l'ownership on free car on rend à l'utilisateur l'interface i_domain / i_part */
  if(part_ext->border_cell_face_idx != NULL) {
    free(part_ext->border_cell_face_idx);
    free(part_ext->border_cell_face);
  }

  if(part_ext->border_face_edge_idx != NULL) {
    free(part_ext->border_face_edge_idx);
    free(part_ext->border_face_edge);
  }

  if(part_ext->border_edge_vtx_idx != NULL) {
    free(part_ext->border_edge_vtx_idx);
    free(part_ext->border_edge_vtx);
  }

  if(part_ext->border_face_vtx_idx != NULL) {
    free(part_ext->border_face_vtx_idx);
    free(part_ext->border_face_vtx);
  }

  if(part_ext->border_face_group_idx != NULL) {
    free(part_ext->border_face_group_idx);
    free(part_ext->border_face_group);
  }

  if(part_ext->border_vtx != NULL) {
    free(part_ext->border_vtx);
  }

  if(part_ext->border_cell_ln_to_gn != NULL) {
    free(part_ext->border_cell_ln_to_gn);
  }

  if(part_ext->border_face_ln_to_gn != NULL) {
    free(part_ext->border_face_ln_to_gn);
  }

  if(part_ext->border_edge_ln_to_gn != NULL) {
    free(part_ext->border_edge_ln_to_gn);
  }

  if(part_ext->border_vtx_ln_to_gn != NULL) {
    free(part_ext->border_vtx_ln_to_gn);
  }

  if(part_ext->border_face_group_ln_to_gn != NULL) {
    free(part_ext->border_face_group_ln_to_gn);
  }

  free(part_ext->n_part_idx);
  free(part_ext->n_part_g_idx);

  for(int i_depth = 0; i_depth < part_ext->depth; ++i_depth) {

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(part_ext->cell_cell_extended_idx         [i_depth][i_part+shift_part]);
        free(part_ext->cell_cell_extended_n           [i_depth][i_part+shift_part]);
        free(part_ext->cell_cell_extended             [i_depth][i_part+shift_part]);
        free(part_ext->unique_order_cell_cell_extended[i_depth][i_part+shift_part]);
        free(part_ext->cell_cell_interface            [i_depth][i_part+shift_part]);
      }
      shift_part += part_ext->n_part[i_domain];
    }
    free(part_ext->cell_cell_extended_idx           [i_depth]);
    free(part_ext->cell_cell_extended_n             [i_depth]);
    free(part_ext->cell_cell_extended               [i_depth]);
    free(part_ext->unique_order_cell_cell_extended  [i_depth]);
    free(part_ext->cell_cell_interface              [i_depth]);
    free(part_ext->n_unique_order_cell_cell_extended[i_depth]);
  }

  free(part_ext->cell_cell_extended_idx);
  free(part_ext->cell_cell_extended_n);
  free(part_ext->cell_cell_extended);
  free(part_ext->unique_order_cell_cell_extended);
  free(part_ext->cell_cell_interface);
  free(part_ext->n_unique_order_cell_cell_extended);

  /* Only shortcut of user data */
  free(part_ext->entity_cell_idx    );
  free(part_ext->entity_cell        );
  free(part_ext->entity_cell_n);

  free(part_ext->dist_neighbor_cell_n   );
  free(part_ext->dist_neighbor_cell_idx );
  free(part_ext->dist_neighbor_cell_desc);
  free(part_ext->dist_neighbor_cell_interface);
  free(part_ext->unique_order_dist_neighbor_cell);
  free(part_ext->n_unique_order_dist_neighbor_cell);

  /* Allocated by distant neightbor */
  free(part_ext->entity_cell_opp_idx);
  free(part_ext->entity_cell_opp_n  );
  free(part_ext->entity_cell_opp    );

  int shift_part = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      free(part_ext->cell_cell_extended_pruned_idx[i_part+shift_part]);
      free(part_ext->cell_cell_extended_pruned    [i_part+shift_part]);
      if(part_ext->owner == PDM_OWNERSHIP_KEEP) {
        free(part_ext->cell_cell_interface_pruned   [i_part+shift_part]);
      }
    }
    shift_part += part_ext->n_part[i_domain];
  }
  free(part_ext->cell_cell_extended_pruned_idx);
  free(part_ext->cell_cell_extended_pruned    );
  free(part_ext->cell_cell_interface_pruned   );

  free(part_ext->n_part);
  part_ext->n_part = NULL;

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(part_ext->parts[i_domain]);
  }
  free(part_ext->parts);

  if(part_ext->shift_by_domain_cell != NULL){
    free(part_ext->shift_by_domain_cell);
  }
  if(part_ext->shift_by_domain_face != NULL){
    free(part_ext->shift_by_domain_face);
  }
  if(part_ext->shift_by_domain_edge != NULL){
    free(part_ext->shift_by_domain_edge);
  }
  if(part_ext->shift_by_domain_vtx != NULL){
    free(part_ext->shift_by_domain_vtx);
  }

  if(part_ext->shift_by_domain_face_group != NULL){
    free(part_ext->shift_by_domain_face_group);
  }

  if(part_ext->owner == PDM_OWNERSHIP_KEEP) {
    if(part_ext->composed_interface_idx != NULL) {
      free(part_ext->composed_interface_idx);
    }

    if(part_ext->composed_interface != NULL) {
      free(part_ext->composed_interface);
    }

    if(part_ext->composed_ln_to_gn_sorted != NULL) {
      free(part_ext->composed_ln_to_gn_sorted);
    }
  }

  free(part_ext);
}


/**
 *
 * \brief Get extended connectivity
 *
 * \param [in]  part_ext            \p PDM_part_extension_t structure instance
 * \param [in]  i_domain            Domain identifier
 * \param [in]  i_part              Partition identifier
 * \param [in]  connectivity_type   Connectivity type
 * \param [out] connect_idx         Connectivity index
 * \param [out] connect             Connectivity
 *
 * \return Number of leading entities
 *
 */

int
PDM_part_extension_connectivity_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                     **connect_idx,
 int                     **connect
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(connectivity_type)
  {
    case PDM_CONNECTIVITY_TYPE_CELL_FACE:
    {
      int n_cell   = part_ext->parts[i_domain][i_part].n_cell;
      n_entity     = part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][n_cell];
      *connect_idx = part_ext->border_cell_face_idx         [shift_part+i_part];
      *connect     = part_ext->border_cell_face             [shift_part+i_part];

      // PDM_log_trace_connectivity_int(*connect_idx, *connect, n_entity, "cell_face");

    }
    break;

    case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
    {
      int n_face   = part_ext->parts[i_domain][i_part].n_face;
      n_entity     = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *connect_idx = part_ext->border_face_edge_idx  [shift_part+i_part];
      *connect     = part_ext->border_face_edge      [shift_part+i_part];
      // PDM_log_trace_connectivity_int(*connect_idx, *connect, n_entity, "face_edge");
    }
    break;

    case PDM_CONNECTIVITY_TYPE_FACE_VTX:
    {
      int n_face   = part_ext->parts[i_domain][i_part].n_face;
      n_entity     = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *connect_idx = part_ext->border_face_vtx_idx   [shift_part+i_part];
      *connect     = part_ext->border_face_vtx       [shift_part+i_part];
    }
    break;

    case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
    {
      int n_edge   = part_ext->parts[i_domain][i_part].n_edge;
      n_entity     = part_ext->edge_edge_extended_idx[shift_part+i_part][n_edge];
      *connect_idx = part_ext->border_edge_vtx_idx   [shift_part+i_part];
      *connect     = part_ext->border_edge_vtx       [shift_part+i_part];
      // PDM_log_trace_connectivity_int(*connect_idx, *connect, n_entity, "edge_vtx");
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

  }

  return n_entity;
}


/**
 *
 * \brief Get global ids
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [out] ln_to_gn     Global ids (size = \ref n_elt)
 *
 * \return  n_elt  Number of elements
 *
 */

int
PDM_part_extension_ln_to_gn_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 PDM_g_num_t             **ln_to_gn
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(mesh_entity)
  {
    case PDM_MESH_ENTITY_CELL:
    {
      int n_cell   = part_ext->parts[i_domain][i_part].n_cell;
      n_entity     = part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][n_cell];
      *ln_to_gn    = part_ext->border_cell_ln_to_gn         [shift_part+i_part];
      // PDM_log_trace_array_int(*ln_to_gn, n_entity, "extent_cell :: ");
    }
    break;

    case PDM_MESH_ENTITY_FACE:
    {
      int n_face   = part_ext->parts[i_domain][i_part].n_face;
      n_entity     = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *ln_to_gn    = part_ext->border_face_ln_to_gn  [shift_part+i_part];
      // PDM_log_trace_array_int(*ln_to_gn, n_entity, "extent_face :: ");
    }
    break;

    case PDM_MESH_ENTITY_EDGE:
    {
      int n_edge   = part_ext->parts[i_domain][i_part].n_edge;
      n_entity     = part_ext->edge_edge_extended_idx[shift_part+i_part][n_edge];
      *ln_to_gn    = part_ext->border_edge_ln_to_gn  [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_VTX:
    {
      int n_vtx   = part_ext->parts[i_domain][i_part].n_vtx;
      n_entity     = part_ext->vtx_vtx_extended_idx[shift_part+i_part][n_vtx];
      *ln_to_gn    = part_ext->border_vtx_ln_to_gn  [shift_part+i_part];
      // PDM_log_trace_array_int(*ln_to_gn, n_entity, "extent_vtx :: ");
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

  }

  return n_entity;
}

/**
 *
 * \brief Get global ids
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [out] ln_to_gn     Global ids (size = \ref n_elt)
 *
 * \return  n_elt  Number of elements
 *
 */

int
PDM_part_extension_interface_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                     **interface_no
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(mesh_entity)
  {
    case PDM_MESH_ENTITY_CELL:
    {
      int n_cell    = part_ext->parts[i_domain][i_part].n_cell;
      n_entity      = part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][n_cell];
      *interface_no = part_ext->cell_cell_interface_pruned   [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_FACE:
    {
      int n_face    = part_ext->parts[i_domain][i_part].n_face;
      n_entity      = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *interface_no = part_ext->face_face_interface   [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_EDGE:
    {
      int n_edge    = part_ext->parts[i_domain][i_part].n_edge;
      n_entity      = part_ext->edge_edge_extended_idx[shift_part+i_part][n_edge];
      *interface_no = part_ext->edge_edge_interface   [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_VTX:
    {
      int n_vtx   =  part_ext->parts[i_domain][i_part].n_vtx;
      n_entity      = part_ext->vtx_vtx_extended_idx[shift_part+i_part][n_vtx];
      *interface_no = part_ext->vtx_vtx_interface   [shift_part+i_part];
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

  }

  return n_entity;
}


int
PDM_part_extension_group_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                     **group_entity_idx,
 int                     **group_entity,
 PDM_g_num_t             **group_entity_ln_to_gn
)
{
  int n_group = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(mesh_entity)
  {
    case PDM_MESH_ENTITY_CELL:
    {
      abort();
    }
    break;

    case PDM_MESH_ENTITY_FACE:
    {
      int n_face_group = part_ext->parts[i_domain][i_part].n_face_group;
      n_group                = n_face_group;
      *group_entity_idx      = part_ext->border_face_group_idx     [shift_part+i_part];
      *group_entity          = part_ext->border_face_group         [shift_part+i_part];
      *group_entity_ln_to_gn = part_ext->border_face_group_ln_to_gn[shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_EDGE:
    {
      abort();
    }
    break;

    case PDM_MESH_ENTITY_VTX:
    {
      abort();
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

  }

  return n_group;
}


/**
 *
 * \brief Get vertex coordinates
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [out] vtx_coord    Vertex coordinates (size = \ref n_vtx * 3)
 *
 * \return  n_vtx  Number of vertices
 *
 */

int
PDM_part_extension_vtx_coord_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 double                  **vtx_coord
)
{
  int shift_part     = part_ext->n_part_idx[i_domain];
  int n_vtx          = part_ext->parts[i_domain][i_part].n_vtx;
  int n_vtx_extended = part_ext->vtx_vtx_extended_idx[shift_part+i_part][n_vtx];
  *vtx_coord = part_ext->border_vtx[shift_part+i_part];

  return n_vtx_extended;
}


int
PDM_part_extension_composed_interface_get
(
 PDM_part_extension_t     *part_ext,
 int                     **composed_interface_idx,
 int                     **composed_interface,
 PDM_g_num_t             **composed_ln_to_gn_sorted
)
{
  *composed_interface_idx   = part_ext->composed_interface_idx;
  *composed_interface       = part_ext->composed_interface;
  *composed_ln_to_gn_sorted = part_ext->composed_ln_to_gn_sorted;

  return part_ext->n_composed_interface;
}


/**
 *
 * \brief Create part_to_part from interior and extended elements
 *
 * \param [out]  ptp                             Part to part structure
 * \param [in]   n_part                          Number of partitions
 * \param [in]   n_int_cell                      Number of interior elements
 * \param [in]   int_cell_ln_to_gn               gnum of interior elements
 * \param [in]   n_ext_cell                      Number of extended elements
 * \param [in]   ext_cell_ln_to_gn               gnum of extended elements
 * \param [out]  n_selected_cell_to_send         Number of elements selected for send
 * \param [out]  selected_cell_to_send           Local numbering of elements selected for send
 *
 */

/* TODO : to generalize for SONICS for vertices
 */

void
PDM_part_to_part_create_from_extension
(
       PDM_part_to_part_t **ptp,
 const int                  n_part,
       int                 *n_int_cell,
 const PDM_g_num_t        **int_cell_ln_to_gn,
       int                 *n_ext_cell,
 const PDM_g_num_t        **ext_cell_ln_to_gn,
       int                **n_selected_cell_to_send,
       int               ***selected_cell_to_send,
 const PDM_MPI_Comm         comm
)
{

  PDM_g_num_t _max_g_num = -1;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    for(int i = 0; i < n_ext_cell[i_part]; ++i) {
      _max_g_num = PDM_MAX(_max_g_num, ext_cell_ln_to_gn[i_part][i]);
    }
    for(int i = 0; i < n_int_cell[i_part]; ++i) {
      _max_g_num = PDM_MAX(_max_g_num, int_cell_ln_to_gn[i_part][i]);
    }
  }

  PDM_g_num_t max_g_num = 0;
  PDM_MPI_Allreduce(&_max_g_num, &max_g_num, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_g_num_t* distrib_cell = PDM_compute_uniform_entity_distribution(comm, max_g_num);

  PDM_part_to_block_t* ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                   PDM_PART_TO_BLOCK_POST_MERGE,
                                                                   1.,
                                             (PDM_g_num_t **)      ext_cell_ln_to_gn,
                                                                   distrib_cell,
                                                                   n_ext_cell,
                                                                   n_part,
                                                                   comm);

  int  block_n_elt = PDM_part_to_block_n_elt_block_get(ptb);
  int *block_n = (int *) malloc( block_n_elt * sizeof(int));
  for(int i = 0; i < block_n_elt; ++i) {
    block_n[i] = 1;
  }

  PDM_g_num_t* distrib_adapt = PDM_part_to_block_adapt_partial_block_to_block(ptb, &block_n, max_g_num);
  free(distrib_adapt);

  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_cell,
                              (const PDM_g_num_t **)  int_cell_ln_to_gn,
                                                      n_int_cell,
                                                      n_part,
                                                      comm);

  int   stride_one = 1;
  int **is_ext_cell = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         block_n,
                         NULL,
           (void ***)    &is_ext_cell);

  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb);
  free(block_n);
  free(distrib_cell);

  int          *_n_selected_cell_to_send        = (int         *)  malloc(n_part * sizeof(int         ));
  int         **_selected_cell_to_send_idx      = (int         **) malloc(n_part * sizeof(int         *));
  int         **_selected_cell_to_send          = (int         **) malloc(n_part * sizeof(int         *));
  PDM_g_num_t **_selected_cell_to_send_ln_to_gn = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    // Compute buffer size
    int n_ext_to_send = 0;
    for(int i = 0; i < n_int_cell[i_part]; ++i) {
      n_ext_to_send += is_ext_cell[i_part][i];
    }

    _n_selected_cell_to_send       [i_part] = n_ext_to_send;
    _selected_cell_to_send_idx     [i_part] = malloc((n_ext_to_send+1) * sizeof(int        ));
    _selected_cell_to_send         [i_part] = malloc( n_ext_to_send    * sizeof(int        ));
    _selected_cell_to_send_ln_to_gn[i_part] = malloc( n_ext_to_send    * sizeof(PDM_g_num_t));

    int idx_write = 0;
    _selected_cell_to_send_idx[i_part][0] = 0;
    for(int i = 0; i < n_int_cell[i_part]; ++i) {
      if(is_ext_cell[i_part][i] == 1) {
        _selected_cell_to_send_idx     [i_part][idx_write+1] = _selected_cell_to_send_idx[i_part][idx_write] + is_ext_cell[i_part][i];
        _selected_cell_to_send         [i_part][idx_write]   = i+1;
        _selected_cell_to_send_ln_to_gn[i_part][idx_write]   = int_cell_ln_to_gn[i_part][i];
        idx_write++;
      }
    }

  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(is_ext_cell[i_part]);
  }
  free(is_ext_cell);

  PDM_part_to_part_t  *_ptp = PDM_part_to_part_create((const PDM_g_num_t **) _selected_cell_to_send_ln_to_gn,
                                                                             _n_selected_cell_to_send,
                                                                             n_part,
                                                      (const PDM_g_num_t **) ext_cell_ln_to_gn,
                                                                             n_ext_cell,
                                                                             n_part,
                                                              (const int **) _selected_cell_to_send_idx,
                                                      (const PDM_g_num_t **) _selected_cell_to_send_ln_to_gn,
                                                                             comm);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(_selected_cell_to_send_idx     [i_part]);
    free(_selected_cell_to_send_ln_to_gn[i_part]);
  }
  free(_selected_cell_to_send_idx);
  free(_selected_cell_to_send_ln_to_gn);

  *n_selected_cell_to_send = _n_selected_cell_to_send;
  *selected_cell_to_send   = _selected_cell_to_send;
  *ptp = _ptp;
}




/**
 *
 * \brief Set connectivity
 *
 * \param [in]  part_ext           \p PDM_part_extension_t structure instance
 * \param [in]  i_domain           Domain identifier
 * \param [in]  i_part             Partition identifier
 * \param [in]  connectivity_type  Type of connectivity
 * \param [in]  connect_idx        Index for connectivity (can be \p NULL for \p PDM_CONNECTIVITY_TYPE_EDGE_VTX)
 * \param [in]  connect            Connectivity
 *
 */

void
PDM_part_extension_connectivity_set
(
 PDM_part_extension_t    *part_ext,
 int                      i_domain,
 int                      i_part,
 PDM_connectivity_type_t  connectivity_type,
 int                     *connect_idx,
 int                     *connect
 )
{
  part_ext->has_connectivity[connectivity_type] = PDM_TRUE;

  switch (connectivity_type) {
    case PDM_CONNECTIVITY_TYPE_EDGE_VTX: {
      part_ext->parts[i_domain][i_part].edge_vtx = connect;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_FACE_VTX: {
      part_ext->parts[i_domain][i_part].face_vtx_idx = connect_idx;
      part_ext->parts[i_domain][i_part].face_vtx     = connect;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_FACE_EDGE: {
      part_ext->parts[i_domain][i_part].face_edge_idx = connect_idx;
      part_ext->parts[i_domain][i_part].face_edge     = connect;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_CELL_FACE: {
      part_ext->parts[i_domain][i_part].cell_face_idx = connect_idx;
      part_ext->parts[i_domain][i_part].cell_face     = connect;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_FACE_CELL: {
      part_ext->parts[i_domain][i_part].face_cell = connect;
      break;
    }

    default: {
      PDM_error(__FILE__, __LINE__, 0, "Connectivity type %d not yet supported\n",
                connectivity_type);
      break;
    }

  }
}


/**
 *
 * \brief Set global ids
 *
 * \param [in]  part_ext     \p PDM_part_extension_t structure instance
 * \param [in]  i_domain     Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [in]  n_entity     Local number of entities
 * \param [in]  ln_to_gn     Global ids (size = \p n_entity)
 *
 */

void
PDM_part_extension_ln_to_gn_set
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                       n_entity,
 PDM_g_num_t              *ln_to_gn
)
{
  switch (mesh_entity) {
    case PDM_MESH_ENTITY_VTX: {
      part_ext->parts[i_domain][i_part].n_vtx         = n_entity;
      part_ext->parts[i_domain][i_part].vtx_ln_to_gn  = ln_to_gn;
      break;
    }

    case PDM_MESH_ENTITY_EDGE: {
      part_ext->parts[i_domain][i_part].n_edge        = n_entity;
      part_ext->parts[i_domain][i_part].edge_ln_to_gn = ln_to_gn;
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      part_ext->parts[i_domain][i_part].n_face        = n_entity;
      part_ext->parts[i_domain][i_part].face_ln_to_gn = ln_to_gn;
      break;
    }

    case PDM_MESH_ENTITY_CELL: {
      part_ext->parts[i_domain][i_part].n_cell        = n_entity;
      part_ext->parts[i_domain][i_part].cell_ln_to_gn = ln_to_gn;
      break;
    }

    default: {
      PDM_error(__FILE__, __LINE__, 0, "Invalid entity type %d\n", mesh_entity);
      break;
    }

  }
}


/**
 *
 * \brief Set vertex coordinates
 *
 * \param [in]  part_ext     \p PDM_part_extension_t structure instance
 * \param [in]  i_domain     Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  vtx_coord    Vertex coordinates (size = 3 * *n_vtx*)
 *
 */

void
PDM_part_extension_vtx_coord_set
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 double                   *vtx_coord
)
{
  part_ext->parts[i_domain][i_part].vtx = vtx_coord;
}


/**
 *
 * \brief Set the connection graph between partitions for the requested entity type
 *
 * \param [in]  multipart             \p PDM_part_extension_t structure instance
 * \param [in]  i_domain              Domain identifier
 * \param [in]  i_part                Partition identifier
 * \param [in]  entity_type           Type of mesh entity
 * \param [in]  part_bound_proc_idx   Partitioning boundary entities index from process (size = *n_rank* + 1)
 * \param [in]  part_bound_part_idx   Partitioning boundary entities index from partition (size = *n_total_part* + 1)
 * \param [in]  part_bound            Partitioning boundary entities (size = 4 * *n_entity_part_bound* = \p part_bound_proc_idx[*n_rank])
 */

void
PDM_part_extension_part_bound_graph_set
(
 PDM_part_extension_t *part_ext,
 int                   i_domain,
 int                   i_part,
 PDM_mesh_entities_t   entity_type,
 int                  *part_bound_proc_idx,
 int                  *part_bound_part_idx,
 int                  *part_bound
)
{
  switch (entity_type) {

    case PDM_MESH_ENTITY_VTX: {
      part_ext->parts[i_domain][i_part].vtx_part_bound_proc_idx  = part_bound_proc_idx;
      part_ext->parts[i_domain][i_part].vtx_part_bound_part_idx  = part_bound_part_idx;
      part_ext->parts[i_domain][i_part].vtx_part_bound           = part_bound;
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      part_ext->parts[i_domain][i_part].face_part_bound_proc_idx = part_bound_proc_idx;
      part_ext->parts[i_domain][i_part].face_part_bound_part_idx = part_bound_part_idx;
      part_ext->parts[i_domain][i_part].face_part_bound          = part_bound;
      break;
    }

    default: {
      PDM_error(__FILE__, __LINE__, 0, "Entity type %d not yet supported\n",
                entity_type);
      break;
    }

  }
}


/**
 *
 * \brief Set group description
 *
 * \param [in]  part_ext               \p PDM_part_extension_t structure instance
 * \param [in]  i_domain               Domain identifier
 * \param [in]  i_part                 Partition identifier
 * \param [in]  entity_type            Type of mesh entity
 * \param [in]  n_group                Number of groups
 * \param [in]  group_entity_idx       Index for group->entity connectivity (size = \p n_group)
 * \param [in]  group_entity           Group->entity connectivity (1-based local ids, size = \p group_entity_idx[\p n_group])
 * \param [in]  group_entity_ln_to_gn  Group->entity connectivity (group-specific global ids, size = \p group_entity_idx[\p n_group])
 *
 */

void
PDM_part_extension_group_set
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                       n_group,
 int                      *group_entity_idx,
 int                      *group_entity,
 PDM_g_num_t              *group_entity_ln_to_gn
)
{
  switch (entity_type) {

    case PDM_MESH_ENTITY_EDGE: {
      part_ext->parts[i_domain][i_part].n_edge_group        = n_group;
      part_ext->parts[i_domain][i_part].edge_bound_idx      = group_entity_idx;
      part_ext->parts[i_domain][i_part].edge_bound          = group_entity;
      part_ext->parts[i_domain][i_part].edge_bound_ln_to_gn = group_entity_ln_to_gn;
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      part_ext->parts[i_domain][i_part].n_face_group        = n_group;
      part_ext->parts[i_domain][i_part].face_bound_idx      = group_entity_idx;
      part_ext->parts[i_domain][i_part].face_bound          = group_entity;
      part_ext->parts[i_domain][i_part].face_bound_ln_to_gn = group_entity_ln_to_gn;
      break;
    }

    default: {
      PDM_error(__FILE__, __LINE__, 0, "Entity type %d not yet supported\n",
                entity_type);
      break;
    }
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

