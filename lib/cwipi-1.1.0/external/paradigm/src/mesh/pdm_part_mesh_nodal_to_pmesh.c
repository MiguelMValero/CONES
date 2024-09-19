
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

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
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_dmesh_nodal_to_dmesh_priv.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_logging.h"
#include "pdm_dmesh.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_unique.h"
#include "pdm_quick_sort.h"
#include "pdm_para_graph_dual.h"
#include "pdm_array.h"
#include "pdm_order.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_mesh_nodal_to_pmesh.h"
#include "pdm_sort.h"
#include "pdm_vtk.h"

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


static
PDM_g_num_t*
_make_absolute_entity_numbering
(
       int          dn_entity,
 const PDM_MPI_Comm comm
)
{
  int n_rank;

  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t n_entity_proc = dn_entity;
  PDM_g_num_t beg_num_abs;

  PDM_MPI_Scan(&n_entity_proc, &beg_num_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  beg_num_abs -= n_entity_proc;

  /** Compute the distribution of elements amont proc **/
  PDM_g_num_t *entity_distrib = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t _dn_face = (PDM_g_num_t) dn_entity;
  PDM_MPI_Allgather((void *) &_dn_face,
                    1,
                    PDM__PDM_MPI_G_NUM,
          (void *) (&entity_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  // entity_distrib[0] = 1;
  entity_distrib[0] = 0;

  for (int i = 1; i < n_rank+1; i++) {
    entity_distrib[i] +=  entity_distrib[i-1];
  }

  return entity_distrib;
}

static
void
_generate_part_entitiy_connectivity
(
  PDM_MPI_Comm   comm,
  int            n_elmt_vol_tot,
  int            n_entity_elt_vol_tot,
  PDM_g_num_t    max_vtx_gnum,
  int           *elmt_entity_vtx_idx,
  PDM_g_num_t   *elmt_entity_vtx,
  int           *elmt_cell_entity_idx,
  int           *parent_elmt_position,
  int           *elmt_entity_kind,
  PDM_g_num_t   *elmt_entity_cell,
  PDM_g_num_t  **dentity_vtx_out,
  int          **dentity_vtx_idx_out,
  PDM_g_num_t  **distrib_entity_out,
  PDM_g_num_t  **elmt_cell_entity_out
)
{
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size (comm, &n_rank);
  PDM_MPI_Comm_rank (comm, &i_rank);

  /*
   * Generate key
   */
  PDM_g_num_t *key_ln_to_gn = malloc(n_entity_elt_vol_tot * sizeof(PDM_g_num_t));
  double      *key_weight   = malloc(n_entity_elt_vol_tot * sizeof(double     ));
  PDM_g_num_t key_mod = 4 * max_vtx_gnum;
  for(int i_entity = 0; i_entity < n_entity_elt_vol_tot; ++i_entity) {
    PDM_g_num_t key = 0;
    for(int idx = elmt_entity_vtx_idx[i_entity]; idx < elmt_entity_vtx_idx[i_entity+1]; ++idx) {
      key += elmt_entity_vtx[idx];
    }
    // min_vtx =
    key_ln_to_gn[i_entity] = key % key_mod + 1;
    key_weight  [i_entity] = elmt_entity_vtx_idx[i_entity+1]-elmt_entity_vtx_idx[i_entity];
  }

  /*
   * Generate distribution
   */
  int sampling_factor = 2;
  int n_iter_max      = 5;
  double tol = 0.10;
  PDM_g_num_t *distrib_key = NULL;
  PDM_distrib_weight(      sampling_factor,
                           n_rank,
                           1,
                           &n_entity_elt_vol_tot,
    (const PDM_g_num_t **) &key_ln_to_gn,
    (const double      **) &key_weight,
                           n_iter_max,
                           tol,
                           comm,
                           &distrib_key);

  free(key_weight);


  /*
   * Prepare send
   */
  int *send_n    = malloc( n_rank              * sizeof(int));
  int *send_idx  = malloc((n_rank+1)           * sizeof(int));
  int *recv_n    = malloc( n_rank              * sizeof(int));
  int *recv_idx  = malloc((n_rank+1)           * sizeof(int));
  int *dest_rank = malloc(n_entity_elt_vol_tot   * sizeof(int));
  for(int i = 0; i < n_rank; ++i) {
    send_n[i] = 0;
  }

  for(int i_entity = 0; i_entity < n_entity_elt_vol_tot; ++i_entity) {
    PDM_g_num_t key = key_ln_to_gn[i_entity];
    int t_rank = PDM_binary_search_gap_long(key-1, distrib_key, n_rank+1);
    dest_rank[i_entity] = t_rank;
    send_n[t_rank]++;
  }

  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, comm);

  int *send_s_entity_vtx_n   = malloc( n_rank    * sizeof(int));
  int *recv_s_entity_vtx_n   = malloc( n_rank    * sizeof(int));
  int *send_s_entity_vtx_idx = malloc((n_rank+1) * sizeof(int));
  int *recv_s_entity_vtx_idx = malloc((n_rank+1) * sizeof(int));
  send_idx[0] = 0;
  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx         [i+1] = send_idx[i] + send_n[i];
    recv_idx         [i+1] = recv_idx[i] + recv_n[i];
    send_n           [i] = 0;
    send_s_entity_vtx_n[i] = 0;
  }

  /*
   * Prepare send
   */
  int         *send_entity_vtx_n     = malloc( n_entity_elt_vol_tot * sizeof(int        ));
  PDM_g_num_t *send_entity_key       = malloc( n_entity_elt_vol_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *send_elmt_entity_cell = malloc( n_entity_elt_vol_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *send_elmt_entity_kind = malloc( n_entity_elt_vol_tot * sizeof(PDM_g_num_t));
  for(int i_entity = 0; i_entity < n_entity_elt_vol_tot; ++i_entity) {
    int t_rank = dest_rank[i_entity];
    int idx_write = send_idx[t_rank] + send_n[t_rank]++;

    send_entity_vtx_n    [idx_write] = elmt_entity_vtx_idx[i_entity+1]-elmt_entity_vtx_idx[i_entity];
    send_entity_key      [idx_write] = key_ln_to_gn  [i_entity];
    send_elmt_entity_cell[idx_write] = elmt_entity_cell[i_entity];
    send_elmt_entity_kind[idx_write] = elmt_entity_kind[i_entity];
    send_s_entity_vtx_n[t_rank] += elmt_entity_vtx_idx[i_entity+1]-elmt_entity_vtx_idx[i_entity];
  }

  free(key_ln_to_gn);
  free(elmt_entity_cell);
  free(elmt_entity_kind);

  int         *recv_entity_vtx_n     = malloc(recv_idx[n_rank] * sizeof(int        ));
  PDM_g_num_t *recv_entity_key       = malloc(recv_idx[n_rank] * sizeof(PDM_g_num_t));
  PDM_g_num_t *recv_elmt_entity_cell = malloc(recv_idx[n_rank] * sizeof(PDM_g_num_t));
  PDM_g_num_t *recv_elmt_entity_kind = malloc(recv_idx[n_rank] * sizeof(PDM_g_num_t));

  PDM_MPI_Alltoallv(send_entity_vtx_n,
                    send_n,
                    send_idx,
                    PDM_MPI_INT,
                    recv_entity_vtx_n,
                    recv_n,
                    recv_idx,
                    PDM_MPI_INT,
                    comm);
  free(send_entity_vtx_n);

  PDM_MPI_Alltoallv(send_entity_key,
                    send_n,
                    send_idx,
                    PDM__PDM_MPI_G_NUM,
                    recv_entity_key,
                    recv_n,
                    recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    comm);
  free(send_entity_key);

  PDM_MPI_Alltoallv(send_elmt_entity_cell,
                    send_n,
                    send_idx,
                    PDM__PDM_MPI_G_NUM,
                    recv_elmt_entity_cell,
                    recv_n,
                    recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    comm);
  free(send_elmt_entity_cell);

  PDM_MPI_Alltoallv(send_elmt_entity_kind,
                    send_n,
                    send_idx,
                    PDM_MPI_INT,
                    recv_elmt_entity_kind,
                    recv_n,
                    recv_idx,
                    PDM_MPI_INT,
                    comm);
  free(send_elmt_entity_kind);

  /*
   * Exchange size of connectivity
   */
  PDM_MPI_Alltoall(send_s_entity_vtx_n, 1, PDM_MPI_INT,
                   recv_s_entity_vtx_n, 1, PDM_MPI_INT, comm);

  send_s_entity_vtx_idx[0] = 0;
  recv_s_entity_vtx_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_s_entity_vtx_idx[i+1] = send_s_entity_vtx_idx[i] + send_s_entity_vtx_n[i];
    recv_s_entity_vtx_idx[i+1] = recv_s_entity_vtx_idx[i] + recv_s_entity_vtx_n[i];
    send_s_entity_vtx_n  [i] = 0;
  }

  /*
   * Exchange of connectivity
   */
  PDM_g_num_t* send_entity_vtx = malloc(send_s_entity_vtx_idx[n_rank] * sizeof(PDM_g_num_t));
  PDM_g_num_t* recv_entity_vtx = malloc(recv_s_entity_vtx_idx[n_rank] * sizeof(PDM_g_num_t));
  for(int i_entity = 0; i_entity < n_entity_elt_vol_tot; ++i_entity) {
    int t_rank = dest_rank[i_entity];
    int idx_write = send_s_entity_vtx_idx[t_rank];

    for(int k = elmt_entity_vtx_idx[i_entity]; k < elmt_entity_vtx_idx[i_entity+1]; ++k) {
      send_entity_vtx[idx_write+send_s_entity_vtx_n[t_rank]++] = elmt_entity_vtx[k];
    }
  }

  PDM_MPI_Alltoallv(send_entity_vtx,
                    send_s_entity_vtx_n,
                    send_s_entity_vtx_idx,
                    PDM__PDM_MPI_G_NUM,
                    recv_entity_vtx,
                    recv_s_entity_vtx_n,
                    recv_s_entity_vtx_idx,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  if(0 == 1) {
    PDM_log_trace_array_int (recv_entity_vtx_n, recv_idx[n_rank]           , "recv_entity_vtx_n ::");
    PDM_log_trace_array_long(recv_entity_key  , recv_idx[n_rank]           , "recv_entity_key   ::");
    PDM_log_trace_array_long(recv_entity_vtx  , recv_s_entity_vtx_idx[n_rank], "recv_entity_vtx   ::");
  }

  int *recv_entity_vtx_idx = malloc((recv_idx[n_rank]+1) * sizeof(int));
  recv_entity_vtx_idx[0] = 0;
  int n_max_vtx = 0;
  for(int i = 0; i < recv_idx[n_rank]; ++i) {
    recv_entity_vtx_idx[i+1] = recv_entity_vtx_idx[i] + recv_entity_vtx_n[i];
    n_max_vtx = PDM_MAX(n_max_vtx, recv_entity_vtx_n[i]);
  }

  /*
   * All data are exchange we need to order the key and resolve conflit in hash table
   */
  int n_recv_key = recv_idx[n_rank];
  int *order = malloc(n_recv_key * sizeof(int));
  PDM_order_gnum_s(recv_entity_key, 1, order, n_recv_key);

  int n_conflit_to_solve = 0;
  PDM_g_num_t last_gnum = -1;

  int *key_conflict_idx = malloc((n_recv_key+1) * sizeof(int));
  key_conflict_idx[0] = 0;
  for(int i = 0; i < n_recv_key; ++i) {
    if(recv_entity_key[order[i]] != last_gnum){
      key_conflict_idx[n_conflit_to_solve+1] = key_conflict_idx[n_conflit_to_solve]+1;
      n_conflit_to_solve++;
      last_gnum = recv_entity_key[order[i]];
    } else {
      key_conflict_idx[n_conflit_to_solve]++;
    }
  }



  int n_max_entity_per_key = 0;
  for(int i = 0; i < n_conflit_to_solve; ++i) {
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, key_conflict_idx[i+1]-key_conflict_idx[i]);
  }


  /*
   * Solve conflict
   */
  if(0 == 1) {
    PDM_log_trace_array_int(key_conflict_idx, n_conflit_to_solve, "key_conflict_idx ::  ");
    for(int i = 0; i < n_conflit_to_solve; ++i) {
      log_trace(" ------ i = %i \n", i);
      for(int i_key = key_conflict_idx[i]; i_key < key_conflict_idx[i+1]; ++i_key) {
        int i_conflict = order[i_key];
        int beg = recv_entity_vtx_idx[i_conflict];
        int n_vtx_in_entity = recv_entity_vtx_idx[i_conflict+1] - beg;
        log_trace(" \t i_key = %i | beg = %i / n = %i : ", recv_entity_key[i_conflict], beg, n_vtx_in_entity);
        PDM_log_trace_array_long(recv_entity_vtx + beg,
                                 n_vtx_in_entity,
                                 "");
      }
    }
    log_trace("-----------------------------------\n\n\n");
  }


  PDM_g_num_t* loc_entity_vtx_1   = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  PDM_g_num_t* loc_entity_vtx_2   = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  int*         already_treat      = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  int*         same_entity_idx    = (int         *) malloc( (n_max_entity_per_key+1) * sizeof(int        ) );
  int*         sens_entity        = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  PDM_g_num_t *tmp_parent         = (PDM_g_num_t *) malloc(3 * n_max_entity_per_key  * sizeof(PDM_g_num_t) );
  int         *order_parent       = (int         *) malloc(n_max_entity_per_key      * sizeof(int        ) );

  /* Resulting array */
  int         *dentity_vtx_idx   = malloc((n_recv_key+1)                * sizeof(int        ));
  PDM_g_num_t *dentity_vtx       = malloc(recv_s_entity_vtx_idx[n_rank] * sizeof(PDM_g_num_t));
  PDM_g_num_t *recv_elmts_entity = malloc( recv_idx[n_rank]             * sizeof(PDM_g_num_t));

  for(int i = 0; i < recv_idx[n_rank] ; ++i) {
    recv_elmts_entity[i] = 0;
  }

  int i_abs_entity   = 0;
  int idx_entity_vtx = 0;
  dentity_vtx_idx[0] = 0;
  for(int i = 0; i < n_conflit_to_solve; ++i) {

    int n_conflict_entitys = key_conflict_idx[i+1] - key_conflict_idx[i];

    for(int j = 0; j < n_conflict_entitys; ++j ) {
      already_treat[j] = -1;
      int i_conflict     = order[key_conflict_idx[i]+j];
      int i_entity1      = order[key_conflict_idx[i]+j];
      int beg_entity_vtx = recv_entity_vtx_idx[i_entity1];
      // int n_vtx_in_entity = recv_entity_vtx_idx[i_entity1+1] - beg_entity_vtx;
      // PDM_g_num_t min_1 = max_vtx_gnum+1;
      // for(int k = 0; k < n_vtx_in_entity; ++k) {
      //   if(recv_entity_vtx[beg_entity_vtx+k] < min_1) {
      //     min_1 = recv_entity_vtx[beg_entity_vtx+k];
      //   };
      // }
      tmp_parent   [3*j  ] = recv_elmt_entity_kind[i_conflict];
      tmp_parent   [3*j+1] = recv_elmt_entity_cell[i_conflict];
      tmp_parent   [3*j+2] = recv_entity_vtx      [beg_entity_vtx];

    }

    PDM_order_gnum_s(tmp_parent, 3, order_parent, n_conflict_entitys);

    if(1 == 0) {
      PDM_log_trace_array_long(tmp_parent  ,2 * n_conflict_entitys, "tmp_parent   ::" );
      PDM_log_trace_array_int (order_parent, n_conflict_entitys, "order_parent ::" );
    }


    for(int idx_entity = 0; idx_entity < n_conflict_entitys; ++idx_entity) {
      int i_entity     = order[key_conflict_idx[i]+order_parent[idx_entity]];

      int beg_1 = recv_entity_vtx_idx[i_entity];
      int n_vtx_in_entity1 = recv_entity_vtx_idx[i_entity+1] - beg_1;

      int idx_next_same_entity = 0;
      sens_entity    [idx_next_same_entity] = 1;
      same_entity_idx[idx_next_same_entity++] = idx_entity;

      /* Only if not treated we search correspondance with other */
      if(already_treat[idx_entity] != 1) {

        PDM_g_num_t key_1 = 0;
        int idx_min_1 = -1;
        PDM_g_num_t min_1 = max_vtx_gnum+1;
        for(int j = 0; j < n_vtx_in_entity1; ++j) {
          loc_entity_vtx_1[j] = recv_entity_vtx[beg_1+j];
          key_1 += loc_entity_vtx_1[j];
          if(loc_entity_vtx_1[j] < min_1) {
            min_1 = loc_entity_vtx_1[j];
            idx_min_1 = j;
          };
        }
        key_1 = key_1 % key_mod + 1;
        PDM_quick_sort_long(loc_entity_vtx_1, 0, n_vtx_in_entity1-1);

        for(int idx_entity2 = 0; idx_entity2 < n_conflict_entitys; ++idx_entity2) {
          int i_entity_next = order[key_conflict_idx[i]+order_parent[idx_entity2]];

          if (i_entity_next == i_entity) {
            continue;
          }

          // printf("conflict : i_entity = %d, i_entity_next = %d...\n", i_entity, i_entity_next);
          if(already_treat[idx_entity2] == 1) {
            continue;
          }

          int beg_2 = recv_entity_vtx_idx[i_entity_next];
          int n_vtx_in_entity2 = recv_entity_vtx_idx[i_entity_next+1] - beg_2;
          if(n_vtx_in_entity1 == n_vtx_in_entity2 ) {

            PDM_g_num_t key_2 = 0;
            int idx_min_2 = -1;
            PDM_g_num_t min_2 = max_vtx_gnum+1;
            for(int j = 0; j < n_vtx_in_entity1; ++j) {
              loc_entity_vtx_2[j] = recv_entity_vtx[beg_2+j];
              key_2 += loc_entity_vtx_2[j];
              if(loc_entity_vtx_2[j] < min_2) {
                min_2 = loc_entity_vtx_2[j];
                idx_min_2 = j;
              };
            }
            key_2 = key_2 % key_mod + 1;
            PDM_quick_sort_long(loc_entity_vtx_2, 0, n_vtx_in_entity2-1);

            // printf("key_1 : %d, key_2 = %d...\n", key_1, key_2);
            if (key_1 != key_2) {
              log_trace("key_1 = %i (beg %d, n %d), key_2 = %i (beg %d)\n",
                        key_1, beg_1, n_vtx_in_entity1, key_2, beg_2);
            }
            assert(key_1 == key_2);

            int is_same_entity = 1;
            for(int i_vtx = 0; i_vtx < n_vtx_in_entity1; ++i_vtx) {
              if(loc_entity_vtx_1[i_vtx] != loc_entity_vtx_2[i_vtx]) {
                is_same_entity = -1;
                break;
              }
            }

            if(is_same_entity == 1 ){
              if(n_vtx_in_entity1 == 2) {
                PDM_g_num_t i1 = recv_entity_vtx[beg_1  ];
                PDM_g_num_t i2 = recv_entity_vtx[beg_1+1];

                PDM_g_num_t j1 = recv_entity_vtx[beg_2];
                PDM_g_num_t j2 = recv_entity_vtx[beg_2+1];
                if(i1 == j1) {
                  sens_entity[idx_next_same_entity] = 1;
                } else {
                  assert(i1 == j2);
                  assert(i2 == j1);
                  sens_entity[idx_next_same_entity] = -1;
                }
              } else {
                // Determine the sens
                PDM_g_num_t i1 = recv_entity_vtx[beg_1 +  idx_min_1                      ];
                PDM_g_num_t i2 = recv_entity_vtx[beg_1 + (idx_min_1+1) % n_vtx_in_entity1];

                PDM_g_num_t j1 = recv_entity_vtx[beg_2 +  idx_min_2                      ];
                PDM_g_num_t j2 = recv_entity_vtx[beg_2 + (idx_min_2+1) % n_vtx_in_entity1];
                // printf(" i1 = %i | i2 = %i | j1 = %i | j2 = %i\n", i1, i2, j1, j2);
                assert(i1 == j1); // Panic
                if(i2 != j2) {
                  sens_entity[idx_next_same_entity] = -1;
                  // printf(" idx3 = %i \n", (idx_min_1+n_vtx_in_entity1) % n_vtx_in_entity1);
                  PDM_g_num_t i3 = recv_entity_vtx[beg_1 + (idx_min_1+n_vtx_in_entity1-1) % n_vtx_in_entity1];
                  // printf(" i1 = %i | i2 = %i | i3 = %i | j1 = %i | j2 = %i\n", i1, i2, i3, j1, j2);
                  assert(i3 == j2);
                } else {
                  sens_entity[idx_next_same_entity] = 1;
                  assert(i2 == j2);
                }
              }
              // Check if same sens
              // same_entity_idx[idx_next_same_entity++] = i_entity_next;
              same_entity_idx[idx_next_same_entity++] = idx_entity2;
            }
          }
        }

        // printf("[%i] - same_entity_idx::", idx_next_same_entity);
        // for(int ii = 0; ii < idx_next_same_entity; ++ii) {
        //   printf(" %i", same_entity_idx[ii]);
        // }
        // printf("\n");

        dentity_vtx_idx[i_abs_entity+1] = dentity_vtx_idx[i_abs_entity] + n_vtx_in_entity1;
        for(int i_vtx = 0; i_vtx < n_vtx_in_entity1; ++i_vtx) {
          dentity_vtx[idx_entity_vtx++] = recv_entity_vtx[beg_1+i_vtx];
        }

        for(int k = 0; k < idx_next_same_entity; ++k) {
          int i_same_entity = same_entity_idx[k];
          int t_entity      = order[key_conflict_idx[i]+order_parent[i_same_entity]];
          int sign = sens_entity[k];
          assert(recv_elmts_entity[t_entity] == 0);
          recv_elmts_entity[t_entity] = sign * (i_abs_entity+1);
          // printf("[%i] - i_same_entity = %i | sign = %i \n", k, i_same_entity, sign);
          already_treat[i_same_entity] = 1;
        }
        i_abs_entity++;

      } /* End if already_treat[i_entity] */
    } /* End conflict */
  }

  if(0 == 1) {
    PDM_log_trace_array_int (order, n_recv_key , "order ::");
    PDM_log_trace_array_long(recv_elmts_entity, n_recv_key , "recv_elmts_entity ::");
    PDM_log_trace_connectivity_long(dentity_vtx_idx, dentity_vtx , i_abs_entity, "dentity_vtx ::");
  }


  free(recv_elmt_entity_cell);
  free(recv_elmt_entity_kind);

  /*
   * Setup distrib
   */
  PDM_g_num_t* entity_distrib = _make_absolute_entity_numbering(i_abs_entity, comm);
  for(int i = 0; i < recv_idx[n_rank]; ++i) {
    PDM_g_num_t g_num = PDM_ABS (recv_elmts_entity[i]);
    int         sgn   = PDM_SIGN(recv_elmts_entity[i]);
    recv_elmts_entity[i] = sgn * (g_num + entity_distrib[i_rank]);
  }


  // PDM_log_trace_array_long(entity_distrib, n_rank+1, "entity_distrib ::");
  /*
   * Reconstruction du entity_ln_to_gn + cell_entity local
   */
  PDM_g_num_t *send_elmt_cell_entity = malloc( elmt_cell_entity_idx[n_elmt_vol_tot] * sizeof(PDM_g_num_t));
  PDM_g_num_t *elmt_cell_entity      = malloc( elmt_cell_entity_idx[n_elmt_vol_tot] * sizeof(PDM_g_num_t));

  PDM_MPI_Alltoallv(recv_elmts_entity,
                    recv_n,
                    recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    send_elmt_cell_entity,
                    send_n,
                    send_idx,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  for(int i = 0; i < n_rank; ++i) {
    send_n[i] = 0;
  }

  for(int i_entity = 0; i_entity < n_entity_elt_vol_tot; ++i_entity) {
    int t_rank = dest_rank[i_entity];
    int idx_write = send_idx[t_rank] + send_n[t_rank]++;
    elmt_cell_entity[i_entity] = send_elmt_cell_entity[idx_write];
  }
  free(send_elmt_cell_entity);

  free(recv_elmts_entity);
  free(loc_entity_vtx_1);
  free(loc_entity_vtx_2);
  free(already_treat   );
  free(same_entity_idx );
  free(sens_entity     );
  free(tmp_parent      );
  free(order_parent    );

  free(key_conflict_idx);
  free(order);

  free(send_n);
  free(send_idx);
  free(recv_n);
  free(recv_idx);
  free(dest_rank);
  free(recv_entity_vtx_n);
  free(recv_entity_vtx_idx);

  free(send_entity_vtx);
  free(recv_entity_vtx);
  free(recv_entity_key);

  free(send_s_entity_vtx_n  );
  free(recv_s_entity_vtx_n  );
  free(send_s_entity_vtx_idx);
  free(recv_s_entity_vtx_idx);

  free(distrib_key);
  free(elmt_entity_vtx_idx     );
  free(elmt_entity_vtx         );
  free(parent_elmt_position  );


  *dentity_vtx_out      = dentity_vtx;
  *dentity_vtx_idx_out  = dentity_vtx_idx;
  *distrib_entity_out   = entity_distrib;
  *elmt_cell_entity_out = elmt_cell_entity;

}

static
PDM_part_mesh_t*
_generate_faces_from_part_mesh_nodal
(
  PDM_part_mesh_nodal_t* pmn
)
{
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size (pmn->comm, &n_rank);
  PDM_MPI_Comm_rank (pmn->comm, &i_rank);

  int n_elmt_vol_tot         = 0;
  int n_face_elt_vol_tot     = 0;
  int n_sum_vtx_vol_face_tot = 0;

  int  n_elmt_surf_tot         = 0;
  int  n_face_elt_surf_tot     = 0;
  int  n_sum_vtx_surf_face_tot = 0;
  int *elmt_face_vtx_idx  = NULL;//malloc((n_face_elt_tot+1) * sizeof(int        ));
  int *elmt_cell_face_idx = NULL;//malloc((n_elmt_tot+1)     * sizeof(int        ));

  int *surf_elmt_face_vtx_idx  = NULL;
  int *surf_elmt_cell_face_idx = NULL;


  PDM_part_mesh_nodal_elmts_decompose_faces_get_size(pmn->volumic,
                                                     &n_elmt_vol_tot,
                                                     &n_face_elt_vol_tot,
                                                     &n_sum_vtx_vol_face_tot,
                                                     &elmt_face_vtx_idx,
                                                     &elmt_cell_face_idx);
  // PDM_log_trace_array_int(elmt_face_vtx_idx, n_face_elt_vol_tot, "elmt_face_vtx_idx : ");
  int have_surface = 0;
  if(pmn->surfacic != NULL) {
    have_surface = 1;
    PDM_part_mesh_nodal_elmts_decompose_faces_get_size(pmn->surfacic,
                                                       &n_elmt_surf_tot,
                                                       &n_face_elt_surf_tot,
                                                       &n_sum_vtx_surf_face_tot,
                                                       &surf_elmt_face_vtx_idx,
                                                       &surf_elmt_cell_face_idx);
    // PDM_log_trace_array_int(surf_elmt_face_vtx_idx, n_face_elt_surf_tot, "surf_elmt_face_vtx_idx : ");
  }


  int n_elmt_tot         = n_elmt_vol_tot         + have_surface * n_elmt_surf_tot;
  int n_face_elt_tot     = n_face_elt_vol_tot     + have_surface * n_face_elt_surf_tot;
  int n_sum_vtx_face_tot = n_sum_vtx_vol_face_tot + have_surface * n_sum_vtx_surf_face_tot;


  PDM_g_num_t **vtx_ln_to_gn = malloc(sizeof(PDM_g_num_t *));
  PDM_g_num_t _max_vtx_gnum = -1;
  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
    vtx_ln_to_gn[i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part);
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
      _max_vtx_gnum = PDM_MAX(_max_vtx_gnum, vtx_ln_to_gn[i_part][i_vtx]);
    }
  }
  PDM_g_num_t max_vtx_gnum = 0;
  PDM_MPI_Allreduce(&_max_vtx_gnum, &max_vtx_gnum, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, pmn->comm);

  // free(elmt_face_vtx_idx );
  // free(elmt_cell_face_idx);
  // elmt_face_vtx_idx  = malloc((n_face_elt_tot+1) * sizeof(int        ));
  // elmt_cell_face_idx = malloc((n_elmt_tot+1)     * sizeof(int        ));

  // int         *elmt_face_vtx_idx    = malloc((n_face_elt_tot+1) * sizeof(int        ));
  PDM_g_num_t *elmt_face_vtx        = malloc(n_sum_vtx_face_tot * sizeof(PDM_g_num_t));
  // int         *elmt_cell_face_idx   = malloc((n_elmt_tot+1)     * sizeof(int        ));
  int         *parent_elmt_position = malloc(n_face_elt_tot     * sizeof(int        ));
  int         *elmt_face_kind       = malloc(n_face_elt_tot     * sizeof(int        ));
  PDM_g_num_t *elmt_face_cell       = malloc(n_face_elt_tot     * sizeof(PDM_g_num_t));

  if(0 == 1) {
    printf("n_face_elt_vol_tot     : %i\n", n_face_elt_vol_tot    );
    printf("n_sum_vtx_vol_face_tot : %i\n", n_sum_vtx_vol_face_tot);
    printf("n_elmt_vol_tot         : %i\n", n_elmt_vol_tot        );

    printf("n_elmt_surf_tot         : %i\n", n_elmt_surf_tot        );
    printf("n_face_elt_surf_tot     : %i\n", n_face_elt_surf_tot    );
    printf("n_sum_vtx_surf_face_tot : %i\n", n_sum_vtx_surf_face_tot);
  }

  // elmt_face_vtx_idx [0] = 0;
  // elmt_cell_face_idx[0] = 0;
  PDM_part_mesh_nodal_elmts_sections_decompose_faces(pmn->volumic,
                                                     vtx_ln_to_gn,
                                                     elmt_face_vtx_idx,
                                                     elmt_face_vtx,
                                                     elmt_cell_face_idx,
                                                     elmt_face_cell,
                                                     parent_elmt_position);


  if (0) {
    // PDM_log_trace_connectivity_long(elmt_face_vtx_idx,
    //                                 elmt_face_vtx,
    //                                 n_face_elt_tot,
    //                                 "elmt_face_vtx : ");

    int *_elmt_face_vtx = malloc(sizeof(int) * n_sum_vtx_face_tot);
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);
    for (int i = 0; i < n_sum_vtx_face_tot; i++) {
      int j = PDM_binary_search_long(elmt_face_vtx[i],
                                     vtx_ln_to_gn[0],
                                     n_vtx);
      assert(j >= 0);
      _elmt_face_vtx[i] = j+1;
    }

    double *coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, 0);

    char filename[999];
    sprintf(filename, "dbg_decomp_faces_%d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           n_vtx,
                           coord,
                           vtx_ln_to_gn[0],
                           n_face_elt_tot,
                           elmt_face_vtx_idx,
                           _elmt_face_vtx,
                           elmt_face_cell,
                           NULL);


    free(_elmt_face_vtx);
  }


  for(int i = 0; i < n_face_elt_vol_tot; ++i) {
    elmt_face_kind[i] = 0;
  }

  if(have_surface == 1) {
    // int         *surf_elmt_face_vtx_idx    = &elmt_face_vtx_idx   [n_face_elt_vol_tot];
    PDM_g_num_t *surf_elmt_face_vtx        = &elmt_face_vtx       [n_sum_vtx_vol_face_tot];
    // int         *surf_elmt_cell_face_idx   = &elmt_cell_face_idx  [n_elmt_vol_tot];
    int         *surf_parent_elmt_position = &parent_elmt_position[n_face_elt_vol_tot];
    int         *surf_elmt_face_kind       = &elmt_face_kind      [n_face_elt_vol_tot];
    PDM_g_num_t *surf_elmt_face_cell       = &elmt_face_cell      [n_face_elt_vol_tot];

    elmt_face_vtx_idx  = realloc(elmt_face_vtx_idx,  (n_face_elt_tot+1) * sizeof(int));
    elmt_cell_face_idx = realloc(elmt_cell_face_idx, (n_elmt_tot+1)     * sizeof(int));
    // int         *surf_elmt_face_vtx_idx   = malloc((n_face_elt_surf_tot+1) * sizeof(int        ));
    // surf_elmt_face_vtx_idx[0] = 0;
    PDM_part_mesh_nodal_elmts_sections_decompose_faces(pmn->surfacic,
                                                       vtx_ln_to_gn,
                                                       surf_elmt_face_vtx_idx,
                                                       surf_elmt_face_vtx,
                                                       surf_elmt_cell_face_idx,
                                                       surf_elmt_face_cell,
                                                       surf_parent_elmt_position);

    for(int i = 0; i < n_elmt_surf_tot; ++i) {
      elmt_cell_face_idx[n_elmt_vol_tot+i+1] = elmt_cell_face_idx[n_elmt_vol_tot+i] + surf_elmt_cell_face_idx[i+1] - surf_elmt_cell_face_idx[i];
    }
    free(surf_elmt_cell_face_idx);

    for(int i = 0; i < n_face_elt_surf_tot; ++i) {
      elmt_face_vtx_idx[n_face_elt_vol_tot+i+1] = elmt_face_vtx_idx[n_face_elt_vol_tot+i] + surf_elmt_face_vtx_idx[i+1] - surf_elmt_face_vtx_idx[i];
    }
    free(surf_elmt_face_vtx_idx);

    for(int i = 0; i < n_face_elt_surf_tot; ++i) {
      surf_elmt_face_kind[i] = 1;
    }

  }

  if(0 == 1) {
    PDM_log_trace_array_int (elmt_face_vtx_idx   , n_face_elt_tot+1  , "elmt_face_vtx_idx    :: ");
    PDM_log_trace_array_int (elmt_cell_face_idx  , n_elmt_tot+1      , "elmt_cell_face_idx   :: ");
    PDM_log_trace_array_int (parent_elmt_position, n_face_elt_tot    , "parent_elmt_position :: ");
    PDM_log_trace_array_int (elmt_face_kind      , n_face_elt_tot    , "elmt_face_kind       :: ");
    PDM_log_trace_array_long(elmt_face_cell      , n_face_elt_tot    , "elmt_face_cell       :: ");
    PDM_log_trace_array_long(elmt_face_vtx       , n_sum_vtx_face_tot, "elmt_face_vtx        :: ");
  }


  PDM_g_num_t *elmt_cell_face  = NULL;
  PDM_g_num_t *entity_distrib  = NULL;
  PDM_g_num_t *dentity_vtx     = NULL;
  int         *dentity_vtx_idx = NULL;
  _generate_part_entitiy_connectivity(pmn->comm,
                                      n_elmt_tot,
                                      n_face_elt_tot,
                                      max_vtx_gnum,
                                      elmt_face_vtx_idx,
                                      elmt_face_vtx,
                                      elmt_cell_face_idx,
                                      parent_elmt_position,
                                      elmt_face_kind,
                                      elmt_face_cell,
                                      &dentity_vtx,
                                      &dentity_vtx_idx,
                                      &entity_distrib,
                                      &elmt_cell_face);


  PDM_part_mesh_t* pm = PDM_part_mesh_create(pmn->n_part,
                                             pmn->comm);

  /*
   * Reconstruction du face_ln_to_gn + cell_face local
   */
  int shift_cell_face = 0;
  int idx_read_elmt   = 0;

  PDM_g_num_t **face_ln_to_gn = malloc(pmn->n_part * sizeof(PDM_g_num_t *));
  int         **cell_face     = malloc(pmn->n_part * sizeof(int         *));
  int          *pn_face       = malloc(pmn->n_part * sizeof(int          ));

  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
    int n_elmts = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_VOLUMIC, i_part);

    int pn_cell_face_idx = elmt_cell_face_idx[idx_read_elmt+n_elmts] - elmt_cell_face_idx[idx_read_elmt];
    face_ln_to_gn[i_part] = malloc(pn_cell_face_idx * sizeof(PDM_g_num_t));
    int *unique_order     = malloc(pn_cell_face_idx * sizeof(int        ));

    int *pcell_face_idx = malloc((n_elmts+1) * sizeof(int));
    pcell_face_idx[0] = 0;
    for(int i = 0; i < n_elmts; ++i) {
      pcell_face_idx[i+1] = pcell_face_idx[i] + elmt_cell_face_idx[idx_read_elmt+i+1] - elmt_cell_face_idx[idx_read_elmt+i];
    }

    PDM_g_num_t *_elmt_cell_face = &elmt_cell_face[shift_cell_face];

    // PDM_log_trace_array_long(_elmt_cell_face, pn_cell_face_idx, "_elmt_cell_face : ");

    for(int i = 0; i < pn_cell_face_idx; ++i) {
      face_ln_to_gn[i_part][i] = PDM_ABS(_elmt_cell_face[i]);
    }

    pn_face      [i_part] = PDM_inplace_unique_long2(face_ln_to_gn[i_part], unique_order, 0, pn_cell_face_idx-1);
    face_ln_to_gn[i_part] = realloc(face_ln_to_gn[i_part], pn_face[i_part] * sizeof(PDM_g_num_t));

    // PDM_log_trace_array_long(face_ln_to_gn[i_part], pn_face[i_part], "face_ln_to_gn ::");

    cell_face[i_part] = malloc(pn_cell_face_idx * sizeof(int));
    for(int idx = 0; idx < pn_cell_face_idx; ++idx) {
      int g_sgn  = PDM_SIGN(_elmt_cell_face[idx]);
      int l_elmt = unique_order[idx];
      cell_face[i_part][idx] = (l_elmt + 1) * g_sgn;
    }

    idx_read_elmt   += n_elmts;
    shift_cell_face += pn_cell_face_idx;

    free(unique_order);

    /*
     * Fill part_mesh
     */
    PDM_part_mesh_n_entity_set(pm,
                               i_part,
                               PDM_MESH_ENTITY_CELL,
                               n_elmts);

    PDM_part_mesh_n_entity_set(pm,
                               i_part,
                               PDM_MESH_ENTITY_FACE,
                               pn_face[i_part]);

    PDM_part_mesh_connectivity_set(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                   cell_face    [i_part],
                                   pcell_face_idx,
                                   PDM_OWNERSHIP_KEEP);


  }

  /*
   * Post-treat surface
   */
  int         **cell_face_bnd     = malloc(pmn->n_part * sizeof(int         *));
  if(have_surface == 1) {
    int n_group = PDM_part_mesh_nodal_n_group_get(pmn, PDM_GEOMETRY_KIND_SURFACIC);
    PDM_part_mesh_n_bound_set(pm, PDM_BOUND_TYPE_FACE, n_group);

    for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
      int n_elmts = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_SURFACIC, i_part);

      int pn_cell_face_bnd_idx = elmt_cell_face_idx[idx_read_elmt+n_elmts] - elmt_cell_face_idx[idx_read_elmt];
      PDM_g_num_t *_elmt_cell_face_bnd = &elmt_cell_face[shift_cell_face];

      cell_face_bnd[i_part] = malloc(pn_cell_face_bnd_idx * sizeof(int));

      for(int idx = 0; idx < pn_cell_face_bnd_idx; ++idx) {
        PDM_g_num_t g_num  = PDM_ABS (_elmt_cell_face_bnd[idx]);
        int g_sgn          = PDM_SIGN(_elmt_cell_face_bnd[idx]);
        int l_elmt         = PDM_binary_search_long(g_num, face_ln_to_gn[i_part], pn_face[i_part]);
        cell_face_bnd[i_part][idx] = (l_elmt + 1) * g_sgn;
      }

      idx_read_elmt   += n_elmts;
      shift_cell_face += pn_cell_face_bnd_idx;

      /*
       * Post-treatment
       */
      for(int i_group = 0; i_group < n_group; ++i_group) {

        int          n_face_group             = 0;
        int         *group_elmt_face          = NULL;
        PDM_g_num_t *group_elmt_face_ln_to_gn = NULL;
        PDM_part_mesh_nodal_group_get(pmn,
                                      PDM_GEOMETRY_KIND_SURFACIC,
                                      i_part,
                                      i_group,
                                      &n_face_group,
                                      &group_elmt_face,
                                      &group_elmt_face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

        int         *group_face          = malloc(n_face_group * sizeof(int        ));
        PDM_g_num_t *group_face_ln_to_gn = malloc(n_face_group * sizeof(PDM_g_num_t));

        for(int idx_face = 0; idx_face < n_face_group; ++idx_face) {
          int i_elmt_face = group_elmt_face[idx_face]-1;
          group_face         [idx_face] = cell_face_bnd[i_part][i_elmt_face];
          group_face_ln_to_gn[idx_face] = group_elmt_face_ln_to_gn[idx_face];
        }

        PDM_part_mesh_bound_set(pm,
                                i_part,
                                i_group,
                                PDM_BOUND_TYPE_FACE,
                                n_face_group,
                                group_face,
                                group_face_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);
      }
    }
  } /* End have_surface == 1 */

  /*
   * Hook face_vtx
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distrib,
                              (const PDM_g_num_t **)  face_ln_to_gn,
                                                      pn_face,
                                                      pmn->n_part,
                                                      pmn->comm);
  int dn_entity = entity_distrib[i_rank+1] - entity_distrib[i_rank];
  int *dentity_vtx_n = malloc(dn_entity * sizeof(int));
  for(int i_entity = 0; i_entity < dn_entity; ++i_entity) {
    dentity_vtx_n[i_entity] = dentity_vtx_idx[i_entity+1] - dentity_vtx_idx[i_entity];
  }

  PDM_g_num_t **pface_vtx_gnum = NULL;
  int         **pface_vtx_n    = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         dentity_vtx_n,
                         dentity_vtx,
                         &pface_vtx_n,
            (void ***)   &pface_vtx_gnum);

  free(dentity_vtx_n);

  int **pface_vtx_idx = malloc(pmn->n_part * sizeof(int *));
  int **pface_vtx     = malloc(pmn->n_part * sizeof(int *));
  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {

    pface_vtx_idx[i_part] = malloc((pn_face[i_part] + 1) * sizeof(int));
    pface_vtx_idx[i_part][0] = 0;
    for(int i_entity = 0; i_entity < pn_face[i_part]; ++i_entity) {
      pface_vtx_idx[i_part][i_entity+1] = pface_vtx_idx[i_part][i_entity] + pface_vtx_n[i_part][i_entity];
    }
    pface_vtx[i_part] = malloc(pface_vtx_idx[i_part][pn_face[i_part]] * sizeof(int));

    /*
     * Translate local
     */
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    PDM_g_num_t  *tmp_vtx_ln_to_gn = malloc(n_vtx * sizeof(PDM_g_num_t));
    int *order_vtx = malloc(n_vtx * sizeof(int));
    for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
      tmp_vtx_ln_to_gn[i_vtx] = vtx_ln_to_gn[i_part][i_vtx];
      order_vtx       [i_vtx] = i_vtx;
    }
    PDM_sort_long(tmp_vtx_ln_to_gn, order_vtx, n_vtx);

    for(int i = 0; i < pface_vtx_idx[i_part][pn_face[i_part]]; ++i) {
      int old_vtx = PDM_binary_search_long(pface_vtx_gnum[i_part][i], tmp_vtx_ln_to_gn, n_vtx);
      pface_vtx[i_part][i] = order_vtx[old_vtx]+1;
    }


    free(tmp_vtx_ln_to_gn);
    free(order_vtx);
    free(pface_vtx_n   [i_part]);
    free(pface_vtx_gnum[i_part]);

    PDM_part_mesh_connectivity_set(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                   pface_vtx    [i_part],
                                   pface_vtx_idx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    PDM_part_mesh_entity_ln_to_gn_set(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      face_ln_to_gn[i_part],
                                      PDM_OWNERSHIP_KEEP);
  }

  free(pface_vtx_n);
  free(pface_vtx_gnum);

  PDM_block_to_part_free(btp);

  if(0 == 1) {
    for(int i_part = 0; i_part < pmn->n_part; ++i_part) {

      int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);
      char filename[999];
      sprintf(filename, "out_pmesh_nodal_to_pmesh_%i_%i.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                             n_vtx,
                             vtx_coord,
                             vtx_ln_to_gn [i_part],
                             pn_face[i_part],
                             pface_vtx_idx[i_part],
                             pface_vtx    [i_part],
                             face_ln_to_gn[i_part],
                             NULL);
    }
  }

  if(have_surface == 1) {
    for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
      free(cell_face_bnd    [i_part]);
    }
  }
  free(cell_face    );
  free(face_ln_to_gn);
  free(pn_face);

  free(cell_face_bnd);
  free(pface_vtx    );
  free(pface_vtx_idx);

  free(elmt_cell_face);
  free(elmt_cell_face_idx);

  free(vtx_ln_to_gn);
  free(dentity_vtx_idx  );
  free(dentity_vtx      );
  free(entity_distrib);

  return pm;
}


static
PDM_part_mesh_t*
_generate_edges_from_part_mesh_nodal
(
  PDM_part_mesh_nodal_t* pmn
)
{
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size (pmn->comm, &n_rank);
  PDM_MPI_Comm_rank (pmn->comm, &i_rank);

  int n_elmt_vol_tot         = 0;
  int n_edge_elt_vol_tot     = 0;
  int n_sum_vtx_vol_edge_tot = 0;

  int n_elmt_surf_tot         = 0;
  int n_edge_elt_surf_tot     = 0;
  int n_sum_vtx_surf_edge_tot = 0;

  PDM_part_mesh_nodal_elmts_decompose_edges_get_size(pmn->surfacic,
                                                     &n_elmt_vol_tot,
                                                     &n_edge_elt_vol_tot,
                                                     &n_sum_vtx_vol_edge_tot );
  int have_ridge = 0;
  if(pmn->ridge != NULL) {
    have_ridge = 1;
    PDM_part_mesh_nodal_elmts_decompose_edges_get_size(pmn->ridge,
                                                       &n_elmt_surf_tot,
                                                       &n_edge_elt_surf_tot,
                                                       &n_sum_vtx_surf_edge_tot);
  }

  int n_elmt_tot         = n_elmt_vol_tot         + have_ridge * n_elmt_surf_tot;
  int n_edge_elt_tot     = n_edge_elt_vol_tot     + have_ridge * n_edge_elt_surf_tot;
  int n_sum_vtx_edge_tot = n_sum_vtx_vol_edge_tot + have_ridge * n_sum_vtx_surf_edge_tot;

  PDM_g_num_t **vtx_ln_to_gn = malloc(sizeof(PDM_g_num_t *));
  PDM_g_num_t _max_vtx_gnum = -1;
  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
    vtx_ln_to_gn[i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part);
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
      _max_vtx_gnum = PDM_MAX(_max_vtx_gnum, vtx_ln_to_gn[i_part][i_vtx]);
    }
  }
  PDM_g_num_t max_vtx_gnum = 0;
  PDM_MPI_Allreduce(&_max_vtx_gnum, &max_vtx_gnum, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, pmn->comm);

  int         *elmt_edge_vtx_idx    = malloc((n_edge_elt_tot+1) * sizeof(int        ));
  PDM_g_num_t *elmt_edge_vtx        = malloc(n_sum_vtx_edge_tot * sizeof(PDM_g_num_t));
  int         *elmt_cell_edge_idx   = malloc((n_elmt_tot+1)     * sizeof(int        ));
  int         *parent_elmt_position = malloc(n_edge_elt_tot     * sizeof(int        ));
  int         *elmt_edge_kind       = malloc(n_edge_elt_tot     * sizeof(int        ));
  PDM_g_num_t *elmt_edge_cell       = malloc(n_edge_elt_tot     * sizeof(PDM_g_num_t));

  if(0 == 1) {
    printf("n_edge_elt_vol_tot     : %i\n", n_edge_elt_vol_tot    );
    printf("n_sum_vtx_vol_edge_tot : %i\n", n_sum_vtx_vol_edge_tot);
    printf("n_elmt_vol_tot         : %i\n", n_elmt_vol_tot        );

    printf("n_elmt_surf_tot         : %i\n", n_elmt_surf_tot        );
    printf("n_edge_elt_surf_tot     : %i\n", n_edge_elt_surf_tot    );
    printf("n_sum_vtx_surf_edge_tot : %i\n", n_sum_vtx_surf_edge_tot);
  }

  elmt_edge_vtx_idx [0] = 0;
  elmt_cell_edge_idx[0] = 0;
  PDM_part_mesh_nodal_elmts_sections_decompose_edges(pmn->surfacic,
                                                     vtx_ln_to_gn,
                                                     elmt_edge_vtx_idx,
                                                     elmt_edge_vtx,
                                                     elmt_cell_edge_idx,
                                                     elmt_edge_cell,
                                                     parent_elmt_position);
  for(int i = 0; i < n_edge_elt_vol_tot; ++i) {
    elmt_edge_kind[i] = 0;
  }

  if(have_ridge == 1) {
    // int         *surf_elmt_face_vtx_idx    = &elmt_face_vtx_idx   [n_face_elt_vol_tot];
    PDM_g_num_t *ridge_elmt_edge_vtx        = &elmt_edge_vtx       [n_sum_vtx_vol_edge_tot];
    int         *ridge_elmt_cell_edge_idx   = &elmt_cell_edge_idx  [n_elmt_vol_tot];
    int         *ridge_parent_elmt_position = &parent_elmt_position[n_edge_elt_vol_tot];
    int         *ridge_elmt_edge_kind       = &elmt_edge_kind      [n_edge_elt_vol_tot];
    PDM_g_num_t *ridge_elmt_edge_cell       = &elmt_edge_cell      [n_edge_elt_vol_tot];

    int         *ridge_elmt_edge_vtx_idx   = malloc((n_edge_elt_surf_tot+1) * sizeof(int        ));
    ridge_elmt_edge_vtx_idx[0] = 0;
    PDM_part_mesh_nodal_elmts_sections_decompose_edges(pmn->ridge,
                                                       vtx_ln_to_gn,
                                                       ridge_elmt_edge_vtx_idx,
                                                       ridge_elmt_edge_vtx,
                                                       ridge_elmt_cell_edge_idx,
                                                       ridge_elmt_edge_cell,
                                                       ridge_parent_elmt_position);
    for(int i = 0; i < n_edge_elt_surf_tot; ++i) {
      elmt_edge_vtx_idx[n_edge_elt_vol_tot+i+1] = elmt_edge_vtx_idx[n_edge_elt_vol_tot+i] + ridge_elmt_edge_vtx_idx[i+1] - ridge_elmt_edge_vtx_idx[i];
    }
    free(ridge_elmt_edge_vtx_idx);

    for(int i = 0; i < n_edge_elt_surf_tot; ++i) {
      ridge_elmt_edge_kind[i] = 1;
    }

  }


  PDM_g_num_t *elmt_cell_edge  = NULL;
  PDM_g_num_t *entity_distrib  = NULL;
  PDM_g_num_t *dentity_vtx     = NULL;
  int         *dentity_vtx_idx = NULL;
  _generate_part_entitiy_connectivity(pmn->comm,
                                      n_elmt_tot,
                                      n_edge_elt_tot,
                                      max_vtx_gnum,
                                      elmt_edge_vtx_idx,
                                      elmt_edge_vtx,
                                      elmt_cell_edge_idx,
                                      parent_elmt_position,
                                      elmt_edge_kind,
                                      elmt_edge_cell,
                                      &dentity_vtx,
                                      &dentity_vtx_idx,
                                      &entity_distrib,
                                      &elmt_cell_edge);

  PDM_part_mesh_t* pm = PDM_part_mesh_create(pmn->n_part,
                                             pmn->comm);

  /*
   * Reconstruction du edge_ln_to_gn + cell_edge local
   */
  int shift_cell_edge = 0;
  int idx_read_elmt   = 0;

  PDM_g_num_t **edge_ln_to_gn = malloc(pmn->n_part * sizeof(PDM_g_num_t *));
  int         **cell_edge     = malloc(pmn->n_part * sizeof(int         *));
  int          *pn_edge       = malloc(pmn->n_part * sizeof(int          ));

  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
    int n_elmts = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_SURFACIC, i_part);

    int pn_cell_edge_idx = elmt_cell_edge_idx[idx_read_elmt+n_elmts] - elmt_cell_edge_idx[idx_read_elmt];
    edge_ln_to_gn[i_part] = malloc(pn_cell_edge_idx * sizeof(PDM_g_num_t));
    int *unique_order     = malloc(pn_cell_edge_idx * sizeof(int        ));

    int *pcell_edge_idx = malloc((n_elmts+1) * sizeof(int));
    pcell_edge_idx[0] = 0;
    for(int i = 0; i < n_elmts; ++i) {
      pcell_edge_idx[i+1] = pcell_edge_idx[i] + elmt_cell_edge_idx[idx_read_elmt+i+1] - elmt_cell_edge_idx[idx_read_elmt+i];
    }

    PDM_g_num_t *_elmt_cell_edge = &elmt_cell_edge[shift_cell_edge];
    for(int i = 0; i < pn_cell_edge_idx; ++i) {
      edge_ln_to_gn[i_part][i] = PDM_ABS(_elmt_cell_edge[i]);
    }

    pn_edge[i_part] = PDM_inplace_unique_long2(edge_ln_to_gn[i_part], unique_order, 0, pn_cell_edge_idx-1);
    edge_ln_to_gn[i_part] = realloc(edge_ln_to_gn[i_part], pn_edge[i_part] * sizeof(PDM_g_num_t));

    // PDM_log_trace_array_long(edge_ln_to_gn[i_part], pn_edge[i_part], "edge_ln_to_gn ::");

    cell_edge[i_part] = malloc(pn_cell_edge_idx * sizeof(int));
    for(int idx = 0; idx < pn_cell_edge_idx; ++idx) {
      int g_sgn  = PDM_SIGN(_elmt_cell_edge[idx]);
      int l_elmt = unique_order[idx];
      cell_edge[i_part][idx] = (l_elmt + 1) * g_sgn;
    }

    idx_read_elmt   += n_elmts;
    shift_cell_edge += pn_cell_edge_idx;

    free(unique_order);

    /*
     * Fill part_mesh
     */
    PDM_part_mesh_n_entity_set(pm,
                               i_part,
                               PDM_MESH_ENTITY_FACE,
                               n_elmts);

    PDM_part_mesh_n_entity_set(pm,
                               i_part,
                               PDM_MESH_ENTITY_EDGE,
                               pn_edge[i_part]);

    PDM_part_mesh_connectivity_set(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                   cell_edge    [i_part],
                                   pcell_edge_idx,
                                   PDM_OWNERSHIP_KEEP);
  }


  /*
   * Post-treat surface
   */
  int         **cell_edge_bnd     = malloc(pmn->n_part * sizeof(int         *));
  if(have_ridge == 1) {
    int n_group = PDM_part_mesh_nodal_n_group_get(pmn, PDM_GEOMETRY_KIND_RIDGE);
    PDM_part_mesh_n_bound_set(pm, PDM_BOUND_TYPE_EDGE, n_group);

    for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
      int n_elmts = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_RIDGE, i_part);

      int pn_cell_edge_bnd_idx = elmt_cell_edge_idx[idx_read_elmt+n_elmts] - elmt_cell_edge_idx[idx_read_elmt];
      PDM_g_num_t *_elmt_cell_edge_bnd = &elmt_cell_edge[shift_cell_edge];

      cell_edge_bnd[i_part] = malloc(pn_cell_edge_bnd_idx * sizeof(int));

      for(int idx = 0; idx < pn_cell_edge_bnd_idx; ++idx) {
        PDM_g_num_t g_num = PDM_ABS (_elmt_cell_edge_bnd[idx]);
        int g_sgn         = PDM_SIGN(_elmt_cell_edge_bnd[idx]);
        int l_elmt        = PDM_binary_search_long(g_num, edge_ln_to_gn[i_part], pn_edge[i_part]);
        cell_edge_bnd[i_part][idx] = (l_elmt + 1) * g_sgn;
      }

      idx_read_elmt   += n_elmts;
      shift_cell_edge += pn_cell_edge_bnd_idx;

      /*
       * Post-treatment
       */
      for(int i_group = 0; i_group < n_group; ++i_group) {

        int          n_edge_group             = 0;
        int         *group_elmt_edge          = NULL;
        PDM_g_num_t *group_elmt_edge_ln_to_gn = NULL;
        PDM_part_mesh_nodal_group_get(pmn,
                                      PDM_GEOMETRY_KIND_RIDGE,
                                      i_part,
                                      i_group,
                                      &n_edge_group,
                                      &group_elmt_edge,
                                      &group_elmt_edge_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

        int         *group_edge          = malloc(n_edge_group * sizeof(int        ));
        PDM_g_num_t *group_edge_ln_to_gn = malloc(n_edge_group * sizeof(PDM_g_num_t));

        for(int idx_edge = 0; idx_edge < n_edge_group; ++idx_edge) {
          int i_elmt_edge = group_elmt_edge[idx_edge]-1;
          group_edge         [idx_edge] = cell_edge_bnd[i_part][i_elmt_edge];
          group_edge_ln_to_gn[idx_edge] = group_elmt_edge_ln_to_gn[idx_edge];
        }

        PDM_part_mesh_bound_set(pm,
                                i_part,
                                i_group,
                                PDM_BOUND_TYPE_EDGE,
                                n_edge_group,
                                group_edge,
                                group_edge_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);
      }
      free(cell_edge_bnd[i_part]);
    }
  } /* End have_ridge == 1 */
  free(cell_edge_bnd);

  /*
   * Hook edge_vtx
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distrib,
                              (const PDM_g_num_t **)  edge_ln_to_gn,
                                                      pn_edge,
                                                      pmn->n_part,
                                                      pmn->comm);
  int dn_entity = entity_distrib[i_rank+1] - entity_distrib[i_rank];
  int *dentity_vtx_n = malloc(dn_entity * sizeof(int));
  for(int i_entity = 0; i_entity < dn_entity; ++i_entity) {
    dentity_vtx_n[i_entity] = dentity_vtx_idx[i_entity+1] - dentity_vtx_idx[i_entity];
  }

  PDM_g_num_t **pedge_vtx_gnum = NULL;
  int         **pedge_vtx_n    = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         dentity_vtx_n,
                         dentity_vtx,
                         &pedge_vtx_n,
            (void ***)   &pedge_vtx_gnum);

  free(dentity_vtx_n);

  int **pedge_vtx_idx = malloc(pmn->n_part * sizeof(int *));
  int **pedge_vtx     = malloc(pmn->n_part * sizeof(int *));
  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {

    pedge_vtx_idx[i_part] = malloc((pn_edge[i_part] + 1) * sizeof(int));
    pedge_vtx_idx[i_part][0] = 0;
    for(int i_entity = 0; i_entity < pn_edge[i_part]; ++i_entity) {
      pedge_vtx_idx[i_part][i_entity+1] = pedge_vtx_idx[i_part][i_entity] + pedge_vtx_n[i_part][i_entity];
    }
    pedge_vtx[i_part] = malloc(pedge_vtx_idx[i_part][pn_edge[i_part]] * sizeof(int));

    /*
     * Translate local
     */
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    PDM_g_num_t  *tmp_vtx_ln_to_gn = malloc(n_vtx * sizeof(PDM_g_num_t));
    int *order_vtx = malloc(n_vtx * sizeof(int));
    for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
      tmp_vtx_ln_to_gn[i_vtx] = vtx_ln_to_gn[i_part][i_vtx];
      order_vtx       [i_vtx] = i_vtx;
    }
    PDM_sort_long(tmp_vtx_ln_to_gn, order_vtx, n_vtx);

    for(int i = 0; i < pedge_vtx_idx[i_part][pn_edge[i_part]]; ++i) {
      int old_vtx = PDM_binary_search_long(pedge_vtx_gnum[i_part][i], tmp_vtx_ln_to_gn, n_vtx);
      pedge_vtx[i_part][i] = order_vtx[old_vtx]+1;
    }


    free(tmp_vtx_ln_to_gn);
    free(order_vtx);
    free(pedge_vtx_n   [i_part]);
    free(pedge_vtx_gnum[i_part]);

    PDM_part_mesh_connectivity_set(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   pedge_vtx    [i_part],
                                   NULL,
                                   PDM_OWNERSHIP_KEEP);

    PDM_part_mesh_entity_ln_to_gn_set(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      edge_ln_to_gn[i_part],
                                      PDM_OWNERSHIP_KEEP);

  }

  free(pedge_vtx_n);
  free(pedge_vtx_gnum);

  PDM_block_to_part_free(btp);

  if(0 == 1) {
    for(int i_part = 0; i_part < pmn->n_part; ++i_part) {

      int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);
      char filename[999];
      sprintf(filename, "out_pmesh_nodal_to_pmesh_%i_%i.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements(filename,
                                 n_vtx,
                                 vtx_coord,
                                 vtx_ln_to_gn [i_part],
                                 PDM_MESH_NODAL_BAR2,
                                 pn_edge[i_part],
                                 pedge_vtx    [i_part],
                                 edge_ln_to_gn[i_part],
                                 0,
                                 NULL,
                                 NULL);
    }
  }

  free(cell_edge    );
  free(edge_ln_to_gn);
  free(pn_edge);

  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
    // free(pedge_vtx    [i_part]);
    free(pedge_vtx_idx[i_part]);
  }
  free(pedge_vtx    );
  free(pedge_vtx_idx);

  free(elmt_cell_edge);
  free(elmt_cell_edge_idx);

  free(vtx_ln_to_gn);
  free(dentity_vtx_idx  );
  free(dentity_vtx      );
  free(entity_distrib);

  return pm;

}


static
PDM_part_mesh_t*
_generate_vtx_from_part_mesh_nodal
(
  PDM_part_mesh_nodal_t* pmn
)
{
  PDM_UNUSED(pmn);
  PDM_error (__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_to_part_mesh with PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_VTX not implemented \n");
  return NULL;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


PDM_part_mesh_t*
PDM_part_mesh_nodal_to_part_mesh
(
        PDM_part_mesh_nodal_t                      *pmn,
  const PDM_dmesh_nodal_to_dmesh_transform_t        transform_kind,
  const PDM_dmesh_nodal_to_dmesh_translate_group_t  transform_group_kind
)
{
  PDM_UNUSED(transform_group_kind);
  PDM_part_mesh_t* pm = NULL;

  switch (transform_kind) {

    case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE:
    {
      pm = _generate_faces_from_part_mesh_nodal(pmn);
    }
    break;

    case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE:
    {
      pm = _generate_edges_from_part_mesh_nodal(pmn);
    }
    break;

    case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_VTX:
    {
      pm = _generate_vtx_from_part_mesh_nodal(pmn);
    }
    break;
  }

  return pm;
}




#ifdef __cplusplus
}
#endif /* __cplusplus */
