
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
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_logging.h"
#include "pdm_dmesh.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_quick_sort.h"
#include "pdm_para_graph_dual.h"
#include "pdm_array.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_sort.h"

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
void
end_timer_and_print(const char* msg, PDM_MPI_Comm comm, double t1){

  double t2 = PDM_MPI_Wtime();

  double delta_t = t2 - t1;
  double delta_max;
  double delta_min;

  PDM_MPI_Allreduce (&delta_t,
                     &delta_max,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     comm);

  PDM_MPI_Allreduce (&delta_t,
                     &delta_min,
                     1,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MIN,
                     comm);

  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);
  if(0 && i_rank == 0) {
    printf("[%i] %s : duration min/max -> %12.5e %12.5e \n", n_rank, msg, delta_min, delta_max);
  }
}

/**
 * \def _compute_keys
 */
static
void
_compute_keys
(
const int          n_face_elt_tot,
const int         *delmt_face_vtx_idx,
const PDM_g_num_t *delmt_face_vtx,
      PDM_g_num_t *ln_to_gn,
      PDM_g_num_t  key_mod
)
{
  for(int i_face = 0; i_face < n_face_elt_tot; ++i_face ) {
    PDM_g_num_t key = 0;
    for(int idx = delmt_face_vtx_idx[i_face]; idx < delmt_face_vtx_idx[i_face+1]; ++idx) {
      key += delmt_face_vtx[idx];
    }
    // min_vtx =
    ln_to_gn[i_face] = key % key_mod + 1;
  }
}

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

  if (0 == 1) {
    printf("beg_num_abs::Face : "PDM_FMT_G_NUM" \n", beg_num_abs);
    printf("entity_distrib : "PDM_FMT_G_NUM,  entity_distrib[0]);
    for (int i = 1; i < n_rank+1; i++) {
      printf(" "PDM_FMT_G_NUM, entity_distrib[i]);
    }
    printf("\n");
  }
  return entity_distrib;
}

static
void
_generate_entitiy_connectivity
(
PDM_MPI_Comm   comm,
PDM_g_num_t   *blk_tot_entity_vtx,
int           *blk_entity_vtx_n,
PDM_g_num_t   *blk_elmt_entity_elmt,
int           *blk_part_id,
int           *blk_parent_elmt_position,
int           *blk_n_entity_per_key,
PDM_g_num_t    n_vtx_abs,
PDM_g_num_t    n_g_child,
int            blk_size,
int            blk_entity_elmt_size,
int            blk_tot_entity_vtx_size,
int            blk_entity_vtx_n_size,
int           *dn_entity,
PDM_g_num_t  **entity_distrib,
int          **dentity_vtx_idx,
PDM_g_num_t  **dentity_vtx,
int          **dentity_elmt_idx,
PDM_g_num_t  **dentity_elmt,
int          **dentity_parent_element_position,
int          **dparent_idx,
PDM_g_num_t  **dparent_gnum,
int          **dparent_sign,
PDM_g_num_t  **delmt_child_distrib,
PDM_g_num_t  **distrib_missing_child,
PDM_g_num_t  **dmissing_child_parent_g_num
)
{

  /*
   * Get the max number of vertex of entitys
   */
  int* blk_entity_vtx_idx  = (int        *) malloc( (blk_entity_vtx_n_size+1) * sizeof(int        ));
  int n_max_entity_per_key = 0;
  int n_tot_entity_per_key = 0;
  int n_child_approx = 0;
  int idx_tmp = 0;
  for(int i_entity = 0; i_entity < blk_size; ++i_entity) {
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, blk_n_entity_per_key[i_entity]);
    n_tot_entity_per_key += blk_n_entity_per_key[i_entity];
    for(int i = 0; i < blk_n_entity_per_key[i_entity]; ++i) {
      if(blk_part_id[idx_tmp] == 1) {
        n_child_approx++;
      }
      idx_tmp++;
    }
  }

  int n_max_vtx       = 0;
  blk_entity_vtx_idx[0] = 0;
  for(int i_entity = 0; i_entity < blk_entity_vtx_n_size; ++i_entity) {
    n_max_vtx          = PDM_MAX(n_max_vtx         , blk_entity_vtx_n    [i_entity]);
    blk_entity_vtx_idx[i_entity+1] = blk_entity_vtx_idx[i_entity] + blk_entity_vtx_n[i_entity];
  }

  // PDM_log_trace_array_int(blk_entity_vtx_idx, blk_entity_vtx_n_size, "blk_entity_vtx_idx:: ");
  /*
   * We need to identify each uniques entitys :
   *      - We have multiple packet to treat
   *      - The connectivity can be sorted in place
   *      - Multiple case can occur :
   *           - Alone entity normaly boundary
   *           - Multiple entitys
   *           - Same entity, we remove and replace by the first
   */
  PDM_g_num_t* loc_entity_vtx_1 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  PDM_g_num_t* loc_entity_vtx_2 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  int*         already_treat    = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  int*         same_entity_idx  = (int         *) malloc( (n_max_entity_per_key+1) * sizeof(int        ) );
  int*         sens_entity      = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  PDM_g_num_t *tmp_parent       = (PDM_g_num_t *) malloc(n_max_entity_per_key      * sizeof(PDM_g_num_t) );
  int         *order            = (int         *) malloc(n_max_entity_per_key      * sizeof(int        ) );

  /*
   * Allocate Memory - entity_vtx - entity_elmt
   */
  *dentity_vtx                     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_tot_entity_vtx_size  );
  *dentity_vtx_idx                 = (int         *) malloc( sizeof(int        ) * (blk_entity_vtx_n_size+1 ));
  *dentity_elmt                    = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  n_tot_entity_per_key     );
  *dentity_elmt_idx                = (int         *) malloc( sizeof(int)         * (blk_entity_elmt_size+1)  );
  *dentity_parent_element_position = (int         *) malloc( sizeof(int)         *  n_tot_entity_per_key     );

  PDM_g_num_t *_dentity_vtx                     = *dentity_vtx;
  int         *_dentity_vtx_idx                 = *dentity_vtx_idx;
  PDM_g_num_t *_dentity_elmt                    = *dentity_elmt;
  int         *_dentity_elmt_idx                = *dentity_elmt_idx;
  // PDM_g_num_t *_dparent_gnum                    = *dparent_gnum;
  int         *_dentity_parent_element_position = *dentity_parent_element_position;
  // printf("blk_entity_elmt_size::%i\n", blk_entity_elmt_size);
  // printf("n_tot_entity_per_key::%i\n", n_tot_entity_per_key);
  // printf("n_child_approx::%i\n", n_child_approx);

  PDM_g_num_t *_tmp_parent_gnum     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  n_child_approx           );
  PDM_g_num_t *_tmp_parent_ln_to_gn = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  n_child_approx           );
  int         *_tmp_parent_sign     = (int         *) malloc( sizeof(int        ) *  n_child_approx           );
  PDM_g_num_t *_tmp_missing_parent_gnum = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  n_child_approx           );
  PDM_g_num_t *_tmp_missing_ln_to_gn = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  n_child_approx           );

  /*
   * Init global numbering
   */
  int i_abs_missing = 0;
  int i_abs_child = 0;
  int i_abs_entity = 0;
  _dentity_vtx_idx[0] = 0;

  int idx = 0;
  int idx_entity_vtx = 0;
  _dentity_elmt_idx[0] = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {
    // printf(" --- Number of conflicting keys :: %i \n", blk_n_entity_per_key[i_key]);

    int n_conflict_entitys = blk_n_entity_per_key[i_key];

    for (int i = 0; i < n_conflict_entitys; i++) {
      tmp_parent[i] = blk_elmt_entity_elmt[idx+i];
      order[i] = i;
    }
    PDM_sort_long(tmp_parent, order, n_conflict_entitys);

    /* Reset */
    for(int j = 0; j < n_conflict_entitys; ++j) {
      already_treat[j] = -1;
    }

    /* Loop over all entitys in conflict */
    //for(int i_entity = 0; i_entity < n_conflict_entitys; ++i_entity) {
    //printf("\nkey %d ------------\n", i_key);
    for (int idx_entity = 0; idx_entity < n_conflict_entitys; ++idx_entity) {
      int i_entity = order[idx_entity];
      //printf("Number of vtx on entitys %i :: %i with index [%i] \n", i_entity, blk_entity_vtx_n[idx+i_entity], idx+i_entity);
      //printf("blk_elmt_entity_elmt[idx+%d] = "PDM_FMT_G_NUM"\n", i_entity, blk_elmt_entity_elmt[idx+i_entity]);

      int n_vtx_entity_1 = blk_entity_vtx_n  [idx+i_entity];
      int beg_1          = blk_entity_vtx_idx[idx+i_entity];
      int idx_next_same_entity = 0;
      sens_entity[idx_next_same_entity] = 1;
      same_entity_idx[idx_next_same_entity++] = i_entity;

      if(already_treat[i_entity] != 1) {

        PDM_g_num_t key_1 = 0;
        int idx_min_1 = -1;
        PDM_g_num_t min_1 = n_vtx_abs+1;
        for(int j = 0; j < n_vtx_entity_1; ++j) {
          loc_entity_vtx_1[j] = blk_tot_entity_vtx[beg_1+j];
          key_1 += loc_entity_vtx_1[j];
          if(loc_entity_vtx_1[j] < min_1) {
            min_1 = loc_entity_vtx_1[j];
            idx_min_1 = j;
          };
        }
        PDM_quick_sort_long(loc_entity_vtx_1, 0, n_vtx_entity_1-1);

        //for(int i_entity_next = i_entity+1; i_entity_next < n_conflict_entitys; ++i_entity_next) {
        for(int idx_entity2 = 0; idx_entity2 < n_conflict_entitys; ++idx_entity2) {
          int i_entity_next = order[idx_entity2];

          if (i_entity_next == i_entity) {
            continue;
          }
          //printf("conflict : i_entity = %d, i_entity_next = %d...\n", i_entity, i_entity_next);
          if(already_treat[i_entity_next] == 1) {
            continue;
          }
          //printf("...not already treated\n");

          int n_vtx_entity_2 = blk_entity_vtx_n[idx+i_entity_next];

          if(n_vtx_entity_1 == n_vtx_entity_2 ) {

            int beg_2 = blk_entity_vtx_idx[idx+i_entity_next];
            PDM_g_num_t key_2 = 0;
            int idx_min_2 = -1;
            PDM_g_num_t min_2 = n_vtx_abs+1;
            for(int j = 0; j < n_vtx_entity_1; ++j) {
              loc_entity_vtx_2[j] = blk_tot_entity_vtx[beg_2+j];
              key_2 += loc_entity_vtx_2[j];
              if(loc_entity_vtx_2[j] < min_2) {
                min_2 = loc_entity_vtx_2[j];
                idx_min_2 = j;
              };
            }
            PDM_quick_sort_long(loc_entity_vtx_2, 0, n_vtx_entity_2-1);

            assert(key_1 == key_2);

            int is_same_entity = 1;
            for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              if(loc_entity_vtx_1[i_vtx] != loc_entity_vtx_2[i_vtx]) {
                is_same_entity = -1;
                break;
              }
            }

            if(is_same_entity == 1 ){
              /*printf("idx_min_1 = %i | idx_min_2 = %i \n", idx_min_1, idx_min_2);
               printf(" test1 :: ");
               for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
                 printf(" %i", (int)blk_tot_entity_vtx[beg_1+i_vtx]);
               }
               printf(" \n");
               printf(" test2 :: ");
               for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
                 printf(" %i", (int)blk_tot_entity_vtx[beg_2+i_vtx]);
               }
               printf(" \n");*/

              if(n_vtx_entity_1 == 2) {
                PDM_g_num_t i1 = blk_tot_entity_vtx[beg_1  ];
                PDM_g_num_t i2 = blk_tot_entity_vtx[beg_1+1];

                PDM_g_num_t j1 = blk_tot_entity_vtx[beg_2];
                PDM_g_num_t j2 = blk_tot_entity_vtx[beg_2+1];
                if(i1 == j1) {
                  sens_entity[idx_next_same_entity] = 1;
                } else {
                  assert(i1 == j2);
                  assert(i2 == j1);
                  sens_entity[idx_next_same_entity] = -1;
                }

              } else {
                // Determine the sens
                PDM_g_num_t i1 = blk_tot_entity_vtx[beg_1 +  idx_min_1                    ];
                PDM_g_num_t i2 = blk_tot_entity_vtx[beg_1 + (idx_min_1+1) % n_vtx_entity_1];

                PDM_g_num_t j1 = blk_tot_entity_vtx[beg_2 +  idx_min_2                    ];
                PDM_g_num_t j2 = blk_tot_entity_vtx[beg_2 + (idx_min_2+1) % n_vtx_entity_1];
                // printf(" i1 = %i | i2 = %i | j1 = %i | j2 = %i\n", i1, i2, j1, j2);
                assert(i1 == j1); // Panic
                if(i2 != j2) {
                  sens_entity[idx_next_same_entity] = -1;
                  // printf(" idx3 = %i \n", (idx_min_1+n_vtx_entity_1) % n_vtx_entity_1);
                  PDM_g_num_t i3 = blk_tot_entity_vtx[beg_1 + (idx_min_1+n_vtx_entity_1-1) % n_vtx_entity_1];
                  // printf(" i1 = %i | i2 = %i | i3 = %i | j1 = %i | j2 = %i\n", i1, i2, i3, j1, j2);
                  assert(i3 == j2);
                } else {
                  sens_entity[idx_next_same_entity] = 1;
                  assert(i2 == j2);
                }
              }

              // Check if same sens
              same_entity_idx[idx_next_same_entity++] = i_entity_next;
            }
          } /* End if same number of vertex */
        } /* End loop next entity */

        /*printf("[%i] - same_entity_idx::", idx_next_same_entity);
          for(int ii = 0; ii < idx_next_same_entity; ++ii) {
            printf(" %i", same_entity_idx[ii]);
          }
          printf("\n");*/
        //PDM_log_trace_array_int(same_entity_idx, idx_next_same_entity, "same_entity_idx : ");

        _dentity_vtx_idx[i_abs_entity+1] = _dentity_vtx_idx[i_abs_entity] + n_vtx_entity_1;
        for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
          _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
          // Ecriture à partir du minumun
          // _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
        }

        _dentity_elmt_idx[i_abs_entity+1] = _dentity_elmt_idx[i_abs_entity];
        int find_child = 0;
        for(int i = 0; i < idx_next_same_entity; ++i) {
          int i_same_entity = same_entity_idx[i];
          int sign = sens_entity[i];

          assert(blk_part_id[idx+i_same_entity] <= 1);
          if(blk_part_id[idx+i_same_entity] == 0) {
            int next_idx = _dentity_elmt_idx[i_abs_entity+1]++;
            _dentity_elmt[next_idx] = sign*blk_elmt_entity_elmt[idx+i_same_entity];
            //printf("[%i] _dentity_elmt[%i] = "PDM_FMT_G_NUM" \n", i_abs_entity, next_idx, sign*blk_elmt_entity_elmt[idx+i_same_entity]);
            _dentity_parent_element_position[next_idx] = blk_parent_elmt_position[idx+i_same_entity];
          } else {
            //printf(" Not treated yet : "PDM_FMT_G_NUM"\n", sign*blk_elmt_entity_elmt[idx+i_same_entity]);
            // _dentity_elmt_child[next_child_idx] = sign*blk_elmt_entity_elmt[idx+i_same_entity];
            // _dentity_parent_element_position_child[next_child_idx] = blk_parent_elmt_position[idx+i_same_entity];
            // printf(" _dparent_gnum[%i] = %i \n", i_abs_child, i_abs_entity );
            if (idx_next_same_entity == 1) {
              /* Dans le cas où on a un elt enfant qui n'est pas retrouvé dans la numérotation
                 parente, on garde l'info qu'il est manquant (utile pour le forçage des bords en génération de maillage */
              _tmp_missing_ln_to_gn[i_abs_missing]    = blk_elmt_entity_elmt[idx+i_same_entity];
              _tmp_missing_parent_gnum[i_abs_missing] = i_abs_entity;
              i_abs_missing++;
            }
            else {
              /* On a trouvé un correspondance entre l'elt enfant et un elt parent, et on garde la correspondance + l'orientation relative */
              _tmp_parent_ln_to_gn[i_abs_child] = blk_elmt_entity_elmt[idx+i_same_entity];
              _tmp_parent_sign    [i_abs_child] = sign;
              //printf(" i_abs_child = %i --> sgn = %i \n", i_abs_child, sign);
              _tmp_parent_gnum[i_abs_child] = i_abs_entity; /* We shift after - We can pass multiple times (ex : edges for quads ) */
              find_child = 1;
            }
          }
          already_treat[i_same_entity] = 1;
        }
        i_abs_child += find_child; /* Only increment if at least one elmet in border in found */
        i_abs_entity++;

      } /* End loop already treated */

    } /* End loop entity in conflict */

    idx += n_conflict_entitys;

  }

  // printf(" realloc dentity_elmt : %i --> %i \n", n_tot_entity_per_key, _dentity_elmt_idx[i_abs_entity]);
  *dentity_elmt = realloc(*dentity_elmt, sizeof(PDM_g_num_t) *  _dentity_elmt_idx[i_abs_entity] );
  _dentity_elmt = *dentity_elmt;

  *dentity_parent_element_position = realloc(*dentity_parent_element_position, sizeof(int) * _dentity_elmt_idx[i_abs_entity] );
  _dentity_parent_element_position = *dentity_parent_element_position;

  /*
   * Free all unused structure
   */
  free(loc_entity_vtx_1);
  free(loc_entity_vtx_2);
  free(blk_entity_vtx_idx);
  free(already_treat);
  free(same_entity_idx);
  free(sens_entity);
  free(blk_tot_entity_vtx);
  free(blk_n_entity_per_key);
  free(blk_entity_vtx_n);
  free(blk_elmt_entity_elmt);
  free(blk_parent_elmt_position);
  free(blk_part_id);

  free (tmp_parent);
  free (order);

  /*
   * Fill up structure
   */
  *dn_entity = i_abs_entity;
  int _dn_entity = *dn_entity;

  /*
   * Realloc
   */
  *dentity_vtx_idx = (int *        ) realloc(*dentity_vtx_idx, (_dn_entity + 1) * sizeof(int * ) );
  _dentity_vtx_idx = *dentity_vtx_idx;

  *dentity_vtx     = (PDM_g_num_t *) realloc(*dentity_vtx    , _dentity_vtx_idx[_dn_entity] * sizeof(PDM_g_num_t * ));
  _dentity_vtx     = *dentity_vtx;

  /*
   * Generate absolute numerotation of entitys
   */
  *entity_distrib = _make_absolute_entity_numbering(_dn_entity, comm);
  PDM_g_num_t* _entity_distrib = *entity_distrib;

  /*
   * Shift child
   */
  int i_rank = -1;
  PDM_MPI_Comm_rank(comm, &i_rank);
  for(int i = 0; i < i_abs_child; ++i) {
    _tmp_parent_gnum[i] = _tmp_parent_gnum[i] + _entity_distrib[i_rank] + 1;
  }
  for(int i = 0; i < i_abs_missing; ++i) {
    _tmp_missing_parent_gnum[i] = _tmp_missing_parent_gnum[i] + _entity_distrib[i_rank] + 1;
  }

  /*
   * Exchange in origin absolute numbering
   */
  double* weights = malloc(i_abs_child * sizeof(double));
  for(int i = 0; i < i_abs_child; ++i) {
    weights[i] = 1.;
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &_tmp_parent_ln_to_gn,
                                                      &weights,
                                                      &i_abs_child,
                                                      1,
                                                      comm);
  free(weights);
  //PDM_g_num_t* distrib = PDM_part_to_block_distrib_index_get(ptb);
  //printf(" A GERE BLOCK PARTIEL INSIDE _generate_entitiy_connectivity in pdm_dmesh_nodal_to_dmesh.c \n");

  int n_rank = -1;
  PDM_MPI_Comm_size(comm, &n_rank);

  int* stride_one = (int *) malloc( i_abs_child * sizeof(int));
  for(int i = 0; i < i_abs_child; ++i) {
    stride_one[i] = 1;
  }

  int* blk_strid = NULL;
  int s_block_data = PDM_part_to_block_exch(ptb,
                                            sizeof(PDM_g_num_t),
                                            PDM_STRIDE_VAR_INTERLACED,
                                            -1,
                                            &stride_one,
                                  (void **) &_tmp_parent_gnum,
                                            &blk_strid,
                                  (void **) dparent_gnum);
  PDM_UNUSED(s_block_data);
  // PDM_log_trace_array_long(*dparent_gnum, s_block_data, "dparent_gnum : ");

  free(blk_strid);
  s_block_data = PDM_part_to_block_exch(ptb,
                                        sizeof(int),
                                        PDM_STRIDE_VAR_INTERLACED,
                                        -1,
                                        &stride_one,
                              (void **) &_tmp_parent_sign,
                                        &blk_strid,
                              (void **) dparent_sign);
  // PDM_log_trace_array_int(*dparent_sign, s_block_data, "dparent_sign : ");

  // log_trace("n_g_child = "PDM_FMT_G_NUM"\n", n_g_child);
  *delmt_child_distrib = PDM_part_to_block_adapt_partial_block_to_block (ptb,
                                                                         &blk_strid,
                                                                         n_g_child);

  int dn_elmt_child = (int) ((*delmt_child_distrib)[i_rank+1] - (*delmt_child_distrib)[i_rank]);
  *dparent_idx = PDM_array_new_idx_from_sizes_int (blk_strid,
                                                   dn_elmt_child);
  // PDM_log_trace_array_int(*dparent_idx, dn_elmt_child+1, "dparent_idx : ");


  PDM_part_to_block_free(ptb);
  free(stride_one);
  free(blk_strid);
  free(_tmp_parent_ln_to_gn);
  free(_tmp_parent_sign);
  free(_tmp_parent_gnum);

  /*
   * Exchange in origin absolute numbering
   */
  // PDM_log_trace_array_long(_tmp_missing_ln_to_gn, i_abs_missing, "_tmp_missing_ln_to_gn : ");
  // PDM_log_trace_array_long(_tmp_missing_parent_gnum, i_abs_missing, "_tmp_missing_parent_gnum : ");

  ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                 PDM_PART_TO_BLOCK_POST_MERGE,
                                 1.,
                                 &_tmp_missing_ln_to_gn,
                                 NULL,
                                 &i_abs_missing,
                                 1,
                                 comm);

  stride_one = (int *) malloc( i_abs_missing * sizeof(int));
  for(int i = 0; i < i_abs_missing; ++i) {
    stride_one[i] = 1;
  }


  s_block_data = PDM_part_to_block_exch(ptb,
                                        sizeof(PDM_g_num_t),
                                        PDM_STRIDE_VAR_INTERLACED,
                                        -1,
                                        &stride_one,
                              (void **) &_tmp_missing_parent_gnum,
                                        &blk_strid,
                              (void **) dmissing_child_parent_g_num);

  int dn_missing_ridge = PDM_part_to_block_n_elt_block_get(ptb);
  // PDM_log_trace_array_int(blk_strid, dn_missing_ridge, "blk_strid : ");
  // PDM_log_trace_array_long(*dmissing_child_parent_g_num, s_block_data, "dmissing_child_parent_g_num : ");

  /*PDM_g_num_t *_distrib_missing_child = PDM_part_to_block_distrib_index_get (ptb);
  *distrib_missing_child = malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  memcpy (*distrib_missing_child, _distrib_missing_child, sizeof(PDM_g_num_t) * (n_rank + 1));*/
  *distrib_missing_child = PDM_compute_entity_distribution (comm, dn_missing_ridge);

  PDM_part_to_block_free(ptb);
  free(stride_one);
  free(blk_strid);
  free(_tmp_missing_ln_to_gn);
  free(_tmp_missing_parent_gnum);

  PDM_g_num_t* _dparent_gnum = *dparent_gnum;
  int*         _dparent_sign = *dparent_sign;

  if( 0 == 1 ){
    printf("i_abs_entity::%i | i_abs_child:: %i \n", i_abs_entity+1, i_abs_child+1);
    PDM_log_trace_array_int(_dentity_vtx_idx , i_abs_entity+1                 , "_dentity_vtx_idx:: " );
    PDM_log_trace_array_long(_dentity_vtx    , _dentity_vtx_idx[i_abs_entity] , "_dentity_vtx:: "     );
    PDM_log_trace_array_int(_dentity_elmt_idx, i_abs_entity+1                 , "_dentity_elmt_idx:: ");
    PDM_log_trace_array_long(_dentity_elmt   , _dentity_elmt_idx[i_abs_entity], "_dentity_elmt:: "    );
    PDM_log_trace_array_long(_dparent_gnum   , i_abs_child                    , "_dparent_gnum:: "    );
    PDM_log_trace_array_int (_dparent_sign   , i_abs_child                    , "_dparent_sign:: "    );
  }

  // end_timer_and_print("PDM_generate_entitiy_connectivity", comm, t1);
}


void
PDM_generate_entitiy_connectivity
(
PDM_MPI_Comm   comm,
PDM_g_num_t    n_vtx_abs,
PDM_g_num_t    n_g_child,
int            n_part,
int           *n_entity_elt_tot,
PDM_g_num_t  **delmt_entity,
int          **delmt_entity_vtx_idx,
PDM_g_num_t  **delmt_entity_vtx,
int          **dparent_elmt_position,
int           *dn_entity,
PDM_g_num_t  **entity_distrib,
int          **dentity_vtx_idx,
PDM_g_num_t  **dentity_vtx,
int          **dentity_elmt_idx,
PDM_g_num_t  **dentity_elmt,
int          **dentity_parent_element_position,
int          **dparent_idx,
PDM_g_num_t  **dparent_gnum,
int          **dparent_sign,
PDM_g_num_t  **delmt_child_distrib,
PDM_g_num_t  **distrib_missing_child,
PDM_g_num_t  **dmissing_child_parent_g_num
)
{
  assert(n_part == 2); // On peux généraliser, et on aura le lien entre tous les niveaux

  PDM_gen_gnum_t* gnum_gen = PDM_gnum_create(3, n_part, PDM_FALSE, 1.e-6, comm, PDM_OWNERSHIP_USER);

  PDM_g_num_t **ln_to_gn           = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  int         **part_id            = (int         **) malloc( n_part * sizeof(int         *));
  double      **weight             = (double      **) malloc( n_part * sizeof(double      *));
  int         **stride_one         = (int         **) malloc( n_part * sizeof(int         *));
  int         **delmt_entity_vtx_n = (int         **) malloc( n_part * sizeof(int         *));
  PDM_g_num_t key_mod = 4 * n_vtx_abs;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    ln_to_gn[i_part] = (PDM_g_num_t *) malloc( n_entity_elt_tot[i_part] * sizeof(PDM_g_num_t));
    part_id [i_part] = (int         *) malloc( n_entity_elt_tot[i_part] * sizeof(int        ));
    _compute_keys(n_entity_elt_tot[i_part],
                  delmt_entity_vtx_idx[i_part],
                  delmt_entity_vtx[i_part],
                  ln_to_gn[i_part],
                  key_mod);

    weight            [i_part] = (double *) malloc( n_entity_elt_tot[i_part] * sizeof(double));
    delmt_entity_vtx_n[i_part] = (int    *) malloc( n_entity_elt_tot[i_part] * sizeof(int   ));
    stride_one        [i_part] = (int    *) malloc( n_entity_elt_tot[i_part] * sizeof(int   ));

    for(int i = 0; i < n_entity_elt_tot[i_part]; ++i) {
      part_id           [i_part][i] = i_part;
      weight            [i_part][i] = 1.;
      stride_one        [i_part][i] = 1.;
      delmt_entity_vtx_n[i_part][i] = delmt_entity_vtx_idx[i_part][i+1] - delmt_entity_vtx_idx[i_part][i];
    }
    free(delmt_entity_vtx_idx[i_part]);

    PDM_gnum_set_from_parents(gnum_gen, i_part, n_entity_elt_tot[i_part], ln_to_gn[i_part]);
  }

  /*
   * Compute an extract gnum
   */
  PDM_gnum_compute(gnum_gen);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(ln_to_gn[i_part]);
    ln_to_gn[i_part] = PDM_gnum_get(gnum_gen, i_part);
    // free(PDM_gnum_get(gnum_gen, i_part));
  }
  PDM_gnum_free(gnum_gen);

  /*
   * Setup part_to_block to filter all keys
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      ln_to_gn,
                                                      weight,
                                                      n_entity_elt_tot,
                                                      n_part,
                                                      comm);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(ln_to_gn[i_part]);
    free(weight  [i_part]);
  }
  free(ln_to_gn);
  free(weight);

  /*
   * Exchange data
   */
  int is_async = 1;

  int         *blk_tot_entity_vtx_n    = NULL;
  PDM_g_num_t *blk_tot_entity_vtx      = NULL;
  int          blk_tot_entity_vtx_size = -1;

  int *blk_n_entity_per_key  = NULL;
  int *blk_entity_vtx_n      = NULL;
  int  blk_entity_vtx_n_size = -1;

  int          *blk_elmt_entity_elmt_stri = NULL;
  PDM_g_num_t  *blk_elmt_entity_elmt      = NULL;
  int           blk_entity_elmt_size      = -1;

  int          *blk_part_id               = NULL;
  int          *blk_parent_elmt_position  = NULL;
  if(is_async == 0) {
    blk_tot_entity_vtx_size = PDM_part_to_block_exch(ptb,
                                                     sizeof(PDM_g_num_t),
                                                     PDM_STRIDE_VAR_INTERLACED,
                                                     -1,
                                                     delmt_entity_vtx_n,
                                           (void **) delmt_entity_vtx,
                                                     &blk_tot_entity_vtx_n,
                                           (void **) &blk_tot_entity_vtx);

    blk_entity_vtx_n_size = PDM_part_to_block_exch(ptb,
                                                   sizeof(int),
                                                   PDM_STRIDE_VAR_INTERLACED,
                                                   -1,
                                                   stride_one,
                                        (void **)  delmt_entity_vtx_n,
                                                  &blk_n_entity_per_key,
                                        (void **) &blk_entity_vtx_n);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(delmt_entity_vtx  [i_part]);
      free(delmt_entity_vtx_n[i_part]);
    }
    free(delmt_entity_vtx_n);

    blk_entity_elmt_size = PDM_part_to_block_exch(ptb,
                                                  sizeof(PDM_g_num_t),
                                                  PDM_STRIDE_VAR_INTERLACED,
                                                  -1,
                                                  stride_one,
                                       (void **)  delmt_entity,
                                                 &blk_elmt_entity_elmt_stri,
                                       (void **) &blk_elmt_entity_elmt);

    free(blk_elmt_entity_elmt_stri); // Same as blk_n_entity_per_key
    PDM_part_to_block_exch(        ptb,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           stride_one,
                (void **)  dparent_elmt_position,
                          &blk_elmt_entity_elmt_stri,
                (void **) &blk_parent_elmt_position);

    free(blk_elmt_entity_elmt_stri); // Same as blk_n_entity_per_key
    PDM_part_to_block_exch(        ptb,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           stride_one,
                (void **)  part_id,
                          &blk_elmt_entity_elmt_stri,
                (void **) &blk_part_id);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(stride_one           [i_part]);
      free(part_id              [i_part]);
      free(delmt_entity         [i_part]);
      free(dparent_elmt_position[i_part]);
    }
    free(stride_one);
    free(part_id);

  } else {

    // PDM_mpi_comm_kind_t k_comm = PDM_MPI_COMM_KIND_WIN_RMA;
    PDM_mpi_comm_kind_t k_comm = PDM_MPI_COMM_KIND_COLLECTIVE;

    int request_entity_vtx = -1;
    PDM_part_to_block_iexch(ptb,
                            k_comm,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            delmt_entity_vtx_n,
                  (void **) delmt_entity_vtx,
                            &blk_tot_entity_vtx_n,
                  (void **) &blk_tot_entity_vtx,
                            &request_entity_vtx);


    int request_entity_vtx_n = -1;
    PDM_part_to_block_iexch(ptb,
                            k_comm,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            stride_one,
                 (void **)  delmt_entity_vtx_n,
                           &blk_n_entity_per_key,
                 (void **) &blk_entity_vtx_n,
                           &request_entity_vtx_n);

    int request_entity_elmt_size = -1;
    PDM_part_to_block_iexch(ptb,
                            k_comm,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            stride_one,
                 (void **)  delmt_entity,
                           &blk_elmt_entity_elmt_stri,
                 (void **) &blk_elmt_entity_elmt,
                           &request_entity_elmt_size);


    int* tmp1_blk_elmt_entity_elmt_stri = NULL;
    int request_parent_elmt_position = -1;
    PDM_part_to_block_iexch(ptb,
                            k_comm,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            stride_one,
                 (void **)  dparent_elmt_position,
                           &tmp1_blk_elmt_entity_elmt_stri,
                 (void **) &blk_parent_elmt_position,
                           &request_parent_elmt_position);

    int* tmp2_blk_elmt_entity_elmt_stri = NULL;
    int request_part_id = -1;
    PDM_part_to_block_iexch(ptb,
                            k_comm,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            stride_one,
                 (void **)  part_id,
                           &tmp2_blk_elmt_entity_elmt_stri,
                 (void **) &blk_part_id,
                           &request_part_id);

    // double dt = PDM_MPI_Wtime() - t1;
    // log_trace("PDM_dmesh_nodal_to_dmesh exch dt = %12.5e \n", dt);
    // t1 = PDM_MPI_Wtime();


    blk_entity_elmt_size   = PDM_part_to_block_iexch_wait(ptb, request_entity_elmt_size);

    PDM_part_to_block_iexch_wait(ptb, request_parent_elmt_position);
    free(tmp1_blk_elmt_entity_elmt_stri);

    PDM_part_to_block_iexch_wait(ptb, request_part_id);
    free(tmp2_blk_elmt_entity_elmt_stri);

    blk_tot_entity_vtx_size = PDM_part_to_block_iexch_wait(ptb, request_entity_vtx);
    blk_entity_vtx_n_size   = PDM_part_to_block_iexch_wait(ptb, request_entity_vtx_n);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(delmt_entity_vtx  [i_part]);
      free(delmt_entity_vtx_n[i_part]);
    }
    free(delmt_entity_vtx_n);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(stride_one           [i_part]);
      free(part_id              [i_part]);
      free(delmt_entity         [i_part]);
      free(dparent_elmt_position[i_part]);
    }
    free(stride_one);
    free(part_id);

    // dt = PDM_MPI_Wtime() - t1;
    // log_trace("PDM_dmesh_nodal_to_dmesh wait dt = %12.5e \n", dt);
  }

  /*
   *  Get the size of the current process bloc
   */
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);

  if( 0 == 1 ) {
    PDM_log_trace_array_int(blk_tot_entity_vtx_n, blk_size             , "blk_tot_entity_vtx_n:: ");
    PDM_log_trace_array_long(blk_tot_entity_vtx , blk_tot_entity_vtx_size, "blk_tot_entity_vtx:: "  );

    PDM_log_trace_array_int(blk_n_entity_per_key, blk_size         , "blk_n_entity_per_key:: ");
    PDM_log_trace_array_int(blk_entity_vtx_n    , blk_entity_vtx_n_size, "blk_entity_vtx_n:: ");

    PDM_log_trace_array_long(blk_elmt_entity_elmt, blk_entity_elmt_size, "blk_elmt_entity_elmt:: ");
    PDM_log_trace_array_int(blk_parent_elmt_position, blk_entity_elmt_size, "blk_parent_elmt_position:: ");
    PDM_log_trace_array_int(blk_part_id, blk_entity_elmt_size, "blk_part_id:: ");
  }

  PDM_part_to_block_free(ptb);
  free(blk_elmt_entity_elmt_stri); // Same as blk_n_entity_per_key
  free(blk_tot_entity_vtx_n);

  _generate_entitiy_connectivity(comm,
                                 blk_tot_entity_vtx,
                                 blk_entity_vtx_n,
                                 blk_elmt_entity_elmt,
                                 blk_part_id,
                                 blk_parent_elmt_position,
                                 blk_n_entity_per_key,
                                 n_vtx_abs,
                                 n_g_child,
                                 blk_size,
                                 blk_entity_elmt_size,
                                 blk_tot_entity_vtx_size,
                                 blk_entity_vtx_n_size,
                                 dn_entity,
                                 entity_distrib,
                                 dentity_vtx_idx,
                                 dentity_vtx,
                                 dentity_elmt_idx,
                                 dentity_elmt,
                                 dentity_parent_element_position,
                                 dparent_idx,
                                 dparent_gnum,
                                 dparent_sign,
                                 delmt_child_distrib,
                                 distrib_missing_child,
                                 dmissing_child_parent_g_num);
}
static
void
_generate_faces_from_dmesh_nodal
(
  _pdm_link_dmesh_nodal_to_dmesh_t *link,
  int                               post_treat_result
)
{

  PDM_dmesh_nodal_t* dmesh_nodal = link->dmesh_nodal;
  assert(dmesh_nodal->volumic  != NULL);
  assert(dmesh_nodal->surfacic != NULL);

  int n_face_elt_vol_tot     = 0;
  int n_sum_vtx_vol_face_tot = 0;

  int n_face_elt_surf_tot     = 0;
  int n_sum_vtx_surf_face_tot = 0;

  PDM_dmesh_nodal_elmts_decompose_faces_get_size(dmesh_nodal->volumic , &n_face_elt_vol_tot , &n_sum_vtx_vol_face_tot );
  PDM_dmesh_nodal_elmts_decompose_faces_get_size(dmesh_nodal->surfacic, &n_face_elt_surf_tot, &n_sum_vtx_surf_face_tot);

  // printf("_generate_faces_from_dmesh_nodal -> n_face_elt_vol_tot      = %i\n", n_face_elt_vol_tot);
  // printf("_generate_faces_from_dmesh_nodal -> n_sum_vtx_vol_face_tot  = %i\n", n_sum_vtx_vol_face_tot);
  // printf("_generate_faces_from_dmesh_nodal -> n_face_elt_surf_tot     = %i\n", n_face_elt_surf_tot);
  // printf("_generate_faces_from_dmesh_nodal -> n_sum_vtx_surf_face_tot = %i\n", n_sum_vtx_surf_face_tot);

  /*
   * Decompose surface
   */
  PDM_g_num_t* delmt_vol_face         = (PDM_g_num_t *) malloc(  n_face_elt_vol_tot     * sizeof(PDM_g_num_t));
  int*         dparent_elmt_vol_pos   = (int         *) malloc(  n_face_elt_vol_tot     * sizeof(int        ));
  int*         delmt_vol_face_vtx_idx = (int         *) malloc( (n_face_elt_vol_tot +1) * sizeof(int        ));
  PDM_g_num_t* delmt_vol_face_vtx     = (PDM_g_num_t *) malloc(  n_sum_vtx_vol_face_tot * sizeof(PDM_g_num_t));

  delmt_vol_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmesh_nodal->volumic,
                               delmt_vol_face_vtx_idx,
                               delmt_vol_face_vtx,
                               delmt_vol_face,
                               NULL, NULL,
                               dparent_elmt_vol_pos);

  /*
   * Decompose surf
   */
  PDM_g_num_t* delmt_surf_face         = (PDM_g_num_t *) malloc(  n_face_elt_surf_tot     * sizeof(PDM_g_num_t));
  int*         dparent_elmt_surf_pos   = (int         *) malloc(  n_face_elt_surf_tot     * sizeof(int        ));
  int*         delmt_surf_face_vtx_idx = (int         *) malloc( (n_face_elt_surf_tot +1) * sizeof(int        ));
  PDM_g_num_t* delmt_surf_face_vtx     = (PDM_g_num_t *) malloc(  n_sum_vtx_surf_face_tot * sizeof(PDM_g_num_t));

  delmt_surf_face_vtx_idx[0] = 0;
  PDM_sections_decompose_faces(dmesh_nodal->surfacic,
                               delmt_surf_face_vtx_idx,
                               delmt_surf_face_vtx,
                               delmt_surf_face,
                               NULL, NULL,
                               dparent_elmt_surf_pos);

  /*
   *  Create empty dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     0,
                                     0,
                                     0,
                                     dmesh_nodal->vtx->n_vtx,
                                     dmesh_nodal->comm);

  /* Juste a view */
  // int dn_vtx = PDM_DMesh_nodal_n_vtx_get(dmesh_nodal);
  // double *dvtx_coords = malloc(3 * dn_vtx * sizeof(double));
  // for(int i_vtx = 0; i_vtx < 3 * n_vtx; ++i_vtx) {
  //   dvtx_coords[i_vtx] = dmesh_nodal->vtx->_coords[i_vtx];
  // }
  // PDM_dmesh_vtx_coord_set(dm, dmesh_nodal->vtx->_coords, PDM_OWNERSHIP_KEEP);
  // dm->_dvtx_coord = dmesh_nodal->vtx->_coords;

  PDM_dmesh_vtx_coord_set(dm, dmesh_nodal->vtx->_coords, PDM_OWNERSHIP_USER);

  assert(link->dmesh == NULL);
  link->dmesh = dm;

  /*
   * Tricks : We use the property of part_to_block to make the concatenation
   *      One part for surfacic
   *      Another for ridge
   */
  int n_part = 2;
  int* n_face_elt_tot = malloc( n_part * sizeof(int));
  n_face_elt_tot[0] = n_face_elt_vol_tot;
  n_face_elt_tot[1] = n_face_elt_surf_tot;

  PDM_g_num_t **delmt_face         = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  int         **dparent_elmt_pos   = (int         **) malloc( n_part * sizeof(int         *));
  int         **delmt_face_vtx_idx = (int         **) malloc( n_part * sizeof(int         *));
  PDM_g_num_t **delmt_face_vtx     = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));

  delmt_face        [0] = delmt_vol_face        ;
  dparent_elmt_pos  [0] = dparent_elmt_vol_pos  ;
  delmt_face_vtx_idx[0] = delmt_vol_face_vtx_idx;
  delmt_face_vtx    [0] = delmt_vol_face_vtx    ;

  delmt_face        [1] = delmt_surf_face        ;
  dparent_elmt_pos  [1] = dparent_elmt_surf_pos  ;
  delmt_face_vtx_idx[1] = delmt_surf_face_vtx_idx;
  delmt_face_vtx    [1] = delmt_surf_face_vtx    ;


  /* Memory is deallocated inside */
  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_CELL] = PDM_TRUE;
  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;



  int n_section_child = dmesh_nodal->surfacic->n_section;
  PDM_g_num_t n_g_child = dmesh_nodal->surfacic->section_distribution[n_section_child];
  if(dmesh_nodal->surfacic->dparent_idx != NULL) {
    free(dmesh_nodal->surfacic->dparent_idx);
  }
  if(dmesh_nodal->surfacic->dparent_gnum != NULL) {
    free(dmesh_nodal->surfacic->dparent_gnum);
  }
  if(dmesh_nodal->surfacic->dparent_sign != NULL) {
    free(dmesh_nodal->surfacic->dparent_sign);
  }
  if(dmesh_nodal->surfacic->delmt_child_distrib != NULL) {
    free(dmesh_nodal->surfacic->delmt_child_distrib);
  }
  PDM_generate_entitiy_connectivity(dmesh_nodal->comm,
                                    dmesh_nodal->n_vtx_abs,
                                    n_g_child,
                                    n_part,
                                    n_face_elt_tot,
                                    delmt_face,
                                    delmt_face_vtx_idx,
                                    delmt_face_vtx,
                                    dparent_elmt_pos,
                                    &dm->dn_face,
                                    &dm->face_distrib,
                                    &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                    &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                    &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_CELL],
                                    &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_CELL],
                                    &link->_dface_parent_element_position,
                                    &dmesh_nodal->surfacic->dparent_idx,
                                    &dmesh_nodal->surfacic->dparent_gnum,
                                    &dmesh_nodal->surfacic->dparent_sign,
                                    &dmesh_nodal->surfacic->delmt_child_distrib,
                                    &link->distrib_missing_surface,
                                    &link->dmissing_surface_parent_g_num);
  free(n_face_elt_tot    );
  free(delmt_face        );
  free(dparent_elmt_pos  );
  free(delmt_face_vtx_idx);
  free(delmt_face_vtx    );


  PDM_g_num_t *_dface_vtx          = dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX];
  int         *_dface_vtx_idx      = dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX];

  PDM_g_num_t* _dface_cell_tmp     = dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_CELL];
  int        * _dface_cell_idx_tmp = dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_CELL];

  // Post_treat
  if (link->distrib_missing_surface[dmesh_nodal->n_rank] == 0 && post_treat_result == 1) {

    // Post_treat
    PDM_g_num_t *dface_cell = (PDM_g_num_t *) malloc( 2 * dm->dn_face * sizeof(PDM_g_num_t));
    int         *dflip_face = (int         *) malloc(     dm->dn_face * sizeof(int        ));
    for(int i_face = 0; i_face < dm->dn_face; ++i_face) {
      dflip_face[i_face] = 1;

      int beg = _dface_cell_idx_tmp[i_face];
      int n_connect_cell = _dface_cell_idx_tmp[i_face+1] - beg;
      if(n_connect_cell == 1) {
        dface_cell[2*i_face  ] = PDM_ABS(_dface_cell_tmp[beg]); // Attention on peut être retourner !!!!
        dface_cell[2*i_face+1] = 0;
      } else {
        assert(n_connect_cell == 2);
        dface_cell[2*i_face  ] = PDM_ABS(_dface_cell_tmp[beg  ]);
        dface_cell[2*i_face+1] = PDM_ABS(_dface_cell_tmp[beg+1]);
      }

      // Flip if the face is in the other sens
      if( PDM_SIGN(_dface_cell_tmp[beg]) == -1 ) {
        dflip_face[i_face] = -1;
        int beg_face_vtx = _dface_vtx_idx[i_face  ];
        int end_face_vtx = _dface_vtx_idx[i_face+1];

        int offset = beg_face_vtx + 1;
        int n_vtx_on_face = end_face_vtx - beg_face_vtx - 1;
        int end_index = n_vtx_on_face - 1;
        PDM_g_num_t tmp_swap;
        for(int i_vtx = 0; i_vtx < n_vtx_on_face/2; ++i_vtx) {
          tmp_swap = _dface_vtx[offset+end_index];
          _dface_vtx[offset+end_index] = _dface_vtx[offset+i_vtx];
          _dface_vtx[offset+i_vtx] = tmp_swap;
          end_index--;
        }
      }
    }


    /*
     * Actualize parent_g_num
     */
    int dn_surfacic = dmesh_nodal->surfacic->delmt_child_distrib[dmesh_nodal->i_rank+1] - dmesh_nodal->surfacic->delmt_child_distrib[dmesh_nodal->i_rank];
    PDM_block_to_part_t* btp = PDM_block_to_part_create(dm->face_distrib,
                                 (const PDM_g_num_t **) &dmesh_nodal->surfacic->dparent_gnum,
                                                        &dn_surfacic,
                                                        1,
                                                        dmesh_nodal->comm);

    int** tmp_pface_flip;
    int stride_one = 1;
    PDM_block_to_part_exch(btp,
                            sizeof(int),
                            PDM_STRIDE_CST_INTERLACED,
                            &stride_one,
               (void *  )   dflip_face,
                            NULL,
               (void ***)  &tmp_pface_flip);
    int *pface_flip = tmp_pface_flip[0];

    // PDM_log_trace_array_int(pface_flip, dn_surfacic, "pface_flip::");
    for(int i = 0; i < dn_surfacic; ++i) {
      dmesh_nodal->surfacic->dparent_sign[i] *= pface_flip[i];
    }
    free(pface_flip);
    free(tmp_pface_flip);
    free(dflip_face);
    PDM_block_to_part_free(btp);

    free(_dface_cell_tmp);
    free(_dface_cell_idx_tmp);
    _dface_cell_tmp     = NULL;
    _dface_cell_idx_tmp = NULL;

    dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_CELL] = dface_cell;
    dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_CELL] = NULL;

    if(0 == 1) {
      PDM_log_trace_array_long(dface_cell, 2 * dm->dn_face, "face_cell::");
      PDM_log_trace_array_long(_dface_vtx, _dface_vtx_idx[dm->dn_face], "_dface_vtx::");
      PDM_log_trace_connectivity_long(_dface_vtx_idx, _dface_vtx, dm->dn_face, "_dface_vtx::");
    }
  }

  // After flip because the orientation of faces changes
  // Not necessary because the algorithim keep the parent to find out boundary condition
  if (link->distrib_missing_surface[dmesh_nodal->n_rank] == 0 && post_treat_result == 1) {
    PDM_setup_connectivity_idx(dm->dn_face,
                               2,
                               dm->dconnectivity[PDM_CONNECTIVITY_TYPE_FACE_CELL],
                               &_dface_cell_idx_tmp,
                               &_dface_cell_tmp);
    for(int i = 0; i < dm->dn_face; ++i) {
      int sgn = 1;
      for(int j = _dface_cell_idx_tmp[i]; j < _dface_cell_idx_tmp[i+1]; ++j){
        _dface_cell_tmp[j] = _dface_cell_tmp[j] * sgn;
        sgn *= -1;
      }
    }
  } else {
    assert(_dface_cell_tmp != NULL);
    assert(_dface_cell_idx_tmp != NULL);
  }

  int is_signed = 1;
  assert(dm->cell_distrib == NULL);
  dm->cell_distrib = (PDM_g_num_t * ) malloc( (dmesh_nodal->n_rank + 1 ) * sizeof(PDM_g_num_t));
  dm->cell_distrib[0] = -1;
  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
  PDM_dconnectivity_transpose(dmesh_nodal->comm,
                              dm->face_distrib,
                              dm->cell_distrib,
                              _dface_cell_idx_tmp,
                              _dface_cell_tmp,
                              is_signed,
                              &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE],
                              &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_CELL_FACE]);

  int dn_cell = dm->cell_distrib[dmesh_nodal->i_rank+1] - dm->cell_distrib[dmesh_nodal->i_rank];
  dm->dn_cell  = dn_cell;
  dm->n_g_cell = dm->cell_distrib[dmesh_nodal->n_rank];
  dm->n_g_face = dm->face_distrib[dmesh_nodal->n_rank];

  if(0 == 1) {
    PDM_log_trace_array_long(dm->dconnectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE], dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE][dn_cell], "dcell_face :: ");
  }

  if (link->distrib_missing_surface[dmesh_nodal->n_rank] == 0 && post_treat_result == 1) {
    free(_dface_cell_idx_tmp);
    free(_dface_cell_tmp);
  }

  int compute_edge_from_faces = 1;
  if( compute_edge_from_faces == 1){
    // Count the number of edges
    int n_edge_face_tot = _dface_vtx_idx[dm->dn_face];

    int*         dface_edge_vtx_idx     = (int         *) malloc( ( n_edge_face_tot + 1) * sizeof(int        ));
    PDM_g_num_t* dface_edge             = (PDM_g_num_t *) malloc(     n_edge_face_tot    * sizeof(PDM_g_num_t));
    PDM_g_num_t* dface_edge_vtx         = (PDM_g_num_t *) malloc( 2 * n_edge_face_tot    * sizeof(PDM_g_num_t));
    int*         dparent_dface_edge_pos = (int         *) malloc(     n_edge_face_tot    * sizeof(int        ));
    int idx = 0;
    int i_edge = 0;
    dface_edge_vtx_idx[0] = 0;
    for(int i_face = 0; i_face < dm->dn_face; ++i_face ) {
      int beg        = _dface_vtx_idx[i_face  ];
      int n_vtx_elmt = _dface_vtx_idx[i_face+1] - beg;

      for(int ivtx = 0; ivtx < n_vtx_elmt; ++ivtx ) {
        int inext = (ivtx + 1) % n_vtx_elmt;

        dface_edge_vtx[idx++] = _dface_vtx[beg+ivtx ];
        dface_edge_vtx[idx++] = _dface_vtx[beg+inext];

        dface_edge_vtx_idx[i_edge+1] = dface_edge_vtx_idx[i_edge] + 2;
        dface_edge[i_edge] = (PDM_g_num_t) i_face + dm->face_distrib[dmesh_nodal->i_rank] + 1;
        dparent_dface_edge_pos[i_edge] = ivtx+1; // Because edge number is same as vtx number in this case
        i_edge++;
      }
    }

    /*
     * Decomposition des ridges si il y a !
     */
    int n_edge_elt_ridge_tot     = 0;
    int n_sum_vtx_ridge_edge_tot = 0;
    if(dmesh_nodal->ridge != NULL) {
      PDM_dmesh_nodal_elmts_decompose_edges_get_size(dmesh_nodal->ridge   , &n_edge_elt_ridge_tot, &n_sum_vtx_ridge_edge_tot);
    }

    /*
     * Decompose ridge
     */
    PDM_g_num_t* delmt_ridge_edge         = (PDM_g_num_t *) malloc(  n_edge_elt_ridge_tot     * sizeof(PDM_g_num_t));
    int*         dparent_elmt_ridge_pos   = (int         *) malloc(  n_edge_elt_ridge_tot     * sizeof(int        ));
    int*         delmt_ridge_edge_vtx_idx = (int         *) malloc( (n_edge_elt_ridge_tot +1) * sizeof(int        ));
    PDM_g_num_t* delmt_ridge_edge_vtx     = (PDM_g_num_t *) malloc(  n_sum_vtx_ridge_edge_tot * sizeof(PDM_g_num_t));

    if(dmesh_nodal->ridge != NULL) {
      delmt_ridge_edge_vtx_idx[0] = 0;
      PDM_sections_decompose_edges(dmesh_nodal->ridge,
                                   delmt_ridge_edge_vtx_idx,
                                   delmt_ridge_edge_vtx,
                                   delmt_ridge_edge,
                                   NULL, NULL,
                                   dparent_elmt_ridge_pos);
    }


    if( 0 == 1 ){
      printf("n_edge_face_tot ::%i\n", n_edge_face_tot );
      PDM_log_trace_array_int (dface_edge_vtx_idx, n_edge_face_tot+1              , "dface_edge_vtx_idx:: ");
      PDM_log_trace_array_long(dface_edge_vtx    , dface_edge_vtx_idx[n_edge_face_tot], "dface_edge_vtx:: ");
      PDM_log_trace_array_long(dface_edge        , n_edge_face_tot                , "dface_edge:: ");
    }

    int n_part_edge = 2;
    int* n_edge_elt_tot = malloc( n_part_edge * sizeof(int));
    n_edge_elt_tot[0] = n_edge_face_tot;
    n_edge_elt_tot[1] = n_edge_elt_ridge_tot;

    PDM_g_num_t **delmt_edge         = (PDM_g_num_t **) malloc( n_part_edge * sizeof(PDM_g_num_t *));
    int         **dparent_face_pos   = (int         **) malloc( n_part_edge * sizeof(int         *));
    int         **delmt_edge_vtx_idx = (int         **) malloc( n_part_edge * sizeof(int         *));
    PDM_g_num_t **delmt_edge_vtx     = (PDM_g_num_t **) malloc( n_part_edge * sizeof(PDM_g_num_t *));

    delmt_edge        [0] = dface_edge;
    dparent_face_pos  [0] = dparent_dface_edge_pos;
    delmt_edge_vtx_idx[0] = dface_edge_vtx_idx;
    delmt_edge_vtx    [0] = dface_edge_vtx    ;

    delmt_edge        [1] = delmt_ridge_edge        ;
    dparent_face_pos  [1] = dparent_elmt_ridge_pos  ;
    delmt_edge_vtx_idx[1] = delmt_ridge_edge_vtx_idx;
    delmt_edge_vtx    [1] = delmt_ridge_edge_vtx    ;

    int n_section_child_ridge = 0;
    PDM_g_num_t n_g_child_ridge = 0;
    if(dmesh_nodal->ridge != NULL) {
      n_section_child_ridge = dmesh_nodal->ridge->n_section;
      n_g_child_ridge       = dmesh_nodal->ridge->section_distribution[n_section_child_ridge];

      if(dmesh_nodal->ridge->dparent_idx != NULL) {
        free(dmesh_nodal->ridge->dparent_idx);
      }
      if(dmesh_nodal->ridge->dparent_gnum != NULL) {
        free(dmesh_nodal->ridge->dparent_gnum);
      }
      if(dmesh_nodal->ridge->dparent_sign != NULL) {
        free(dmesh_nodal->ridge->dparent_sign);
      }
      if(dmesh_nodal->ridge->delmt_child_distrib != NULL) {
        free(dmesh_nodal->ridge->delmt_child_distrib);
      }

      PDM_generate_entitiy_connectivity(dmesh_nodal->comm,
                                        dmesh_nodal->n_vtx_abs,
                                        n_g_child_ridge,
                                        n_part_edge,
                                        n_edge_elt_tot,
                                        delmt_edge,
                                        delmt_edge_vtx_idx,
                                        delmt_edge_vtx,
                                        dparent_face_pos,
                                        &dm->dn_edge,
                                        &dm->edge_distrib,
                                        &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                        &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                        &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                        &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                        &link->_dedge_parent_element_position,
                                        &dmesh_nodal->ridge->dparent_idx,
                                        &dmesh_nodal->ridge->dparent_gnum,
                                        &dmesh_nodal->ridge->dparent_sign,
                                        &dmesh_nodal->ridge->delmt_child_distrib,
                                        &link->distrib_missing_ridge,
                                        &link->dmissing_ridge_parent_g_num);
      dm->n_g_edge = dm->edge_distrib[dmesh_nodal->n_rank];
    } else {
      int         *dparent_idx = NULL;
      PDM_g_num_t *dparent_gnum = NULL;
      int         *dparent_sign = NULL;
      PDM_g_num_t *delmt_child_distrib = NULL;
      PDM_generate_entitiy_connectivity(dmesh_nodal->comm,
                                        dmesh_nodal->n_vtx_abs,
                                        n_g_child_ridge,
                                        n_part_edge,
                                        n_edge_elt_tot,
                                        delmt_edge,
                                        delmt_edge_vtx_idx,
                                        delmt_edge_vtx,
                                        dparent_face_pos,
                                        &dm->dn_edge,
                                        &dm->edge_distrib,
                                        &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                        &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                        &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                        &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                        &link->_dedge_parent_element_position,
                                        &dparent_idx,
                                        &dparent_gnum,
                                        &dparent_sign,
                                        &delmt_child_distrib,
                                        &link->distrib_missing_ridge,
                                        &link->dmissing_ridge_parent_g_num);
      dm->n_g_edge = dm->edge_distrib[dmesh_nodal->n_rank];

      free(dparent_idx);
      free(dparent_gnum);
      free(dparent_sign);
      free(delmt_child_distrib);
    }
    free(n_edge_elt_tot    );
    free(delmt_edge        );
    free(dparent_face_pos  );
    free(delmt_edge_vtx_idx);
    free(delmt_edge_vtx    );

    dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_FACE] = PDM_TRUE;

    if( 0 == 1 ){
      printf("dmesh_nodal->dn_edge ::%i\n", dm->dn_edge );
      PDM_log_trace_array_int (dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX], dm->dn_edge+1                  , "dm->_dedge_vtx_idx:: ");
      PDM_log_trace_array_long(dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX], dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX][dm->dn_edge], "dm->_dedge_vtx:: ");

      // PDM_log_trace_array_int (dm->dface_edge_idx, dm->dn_face+1                           , "dm->dface_edge_idx:: ");
      // PDM_log_trace_array_long(dm->dface_edge    , dm->dface_edge_idx[dm->dn_face], "dm->dface_edge:: ");

      PDM_log_trace_array_int (dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE], dm->dn_edge+1                   , "dm->_dedge_face_idx:: ");
      PDM_log_trace_array_long(dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_FACE], dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE][dm->dn_edge], "dm->_dedge_face:: ");
    }

    assert(dm->edge_distrib != NULL);
    assert(dm->face_distrib != NULL);
    dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    PDM_dconnectivity_transpose(dmesh_nodal->comm,
                                dm->edge_distrib,
                                dm->face_distrib,
                                dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                1, // is_signed
                                &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE]);

    // PDM_log_trace_array_int (dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE], dm->dn_face+1                   , "dm->_dface_edge_idx:: ");
    // PDM_log_trace_array_long(dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE], dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][dm->dn_face], "dm->_dface_edge:: ");

  }
}

static
void
_generate_edges_from_dmesh_nodal
(
  _pdm_link_dmesh_nodal_to_dmesh_t *link,
  int                               post_treat_result
)
{
  PDM_dmesh_nodal_t* dmesh_nodal = link->dmesh_nodal;
  assert(dmesh_nodal->surfacic != NULL);
  assert(dmesh_nodal->ridge    != NULL);

  int n_edge_elt_surf_tot     = 0;
  int n_sum_vtx_surf_edge_tot = 0;

  int n_edge_elt_ridge_tot     = 0;
  int n_sum_vtx_ridge_edge_tot = 0;

  PDM_dmesh_nodal_elmts_decompose_edges_get_size(dmesh_nodal->surfacic, &n_edge_elt_surf_tot , &n_sum_vtx_surf_edge_tot );
  PDM_dmesh_nodal_elmts_decompose_edges_get_size(dmesh_nodal->ridge   , &n_edge_elt_ridge_tot, &n_sum_vtx_ridge_edge_tot);

  // printf("_generate_edges_from_dmesh_nodal -> n_edge_elt_surf_tot      = %i\n", n_edge_elt_surf_tot);
  // printf("_generate_edges_from_dmesh_nodal -> n_sum_vtx_surf_edge_tot  = %i\n", n_sum_vtx_surf_edge_tot);
  // printf("_generate_edges_from_dmesh_nodal -> n_edge_elt_ridge_tot     = %i\n", n_edge_elt_ridge_tot);
  // printf("_generate_edges_from_dmesh_nodal -> n_sum_vtx_ridge_edge_tot = %i\n", n_sum_vtx_ridge_edge_tot);

  /*
   * Decompose surface
   */
  PDM_g_num_t* delmt_surf_edge         = (PDM_g_num_t *) malloc(  n_edge_elt_surf_tot     * sizeof(PDM_g_num_t));
  int*         dparent_elmt_surf_pos   = (int         *) malloc(  n_edge_elt_surf_tot     * sizeof(int        ));
  int*         delmt_surf_edge_vtx_idx = (int         *) malloc( (n_edge_elt_surf_tot +1) * sizeof(int        ));
  PDM_g_num_t* delmt_surf_edge_vtx     = (PDM_g_num_t *) malloc(  n_sum_vtx_surf_edge_tot * sizeof(PDM_g_num_t));


  delmt_surf_edge_vtx_idx[0] = 0;
  PDM_sections_decompose_edges(dmesh_nodal->surfacic,
                               delmt_surf_edge_vtx_idx,
                               delmt_surf_edge_vtx,
                               delmt_surf_edge,
                               NULL, NULL,
                               dparent_elmt_surf_pos);

  if (0) {
    PDM_log_trace_connectivity_long(delmt_surf_edge_vtx_idx, delmt_surf_edge_vtx, n_edge_elt_surf_tot, "delmt_surf_edge_vtx : ");
  }

  /*
   * Decompose ridge
   */
  PDM_g_num_t* delmt_ridge_edge         = (PDM_g_num_t *) malloc(  n_edge_elt_ridge_tot     * sizeof(PDM_g_num_t));
  int*         dparent_elmt_ridge_pos   = (int         *) malloc(  n_edge_elt_ridge_tot     * sizeof(int        ));
  int*         delmt_ridge_edge_vtx_idx = (int         *) malloc( (n_edge_elt_ridge_tot +1) * sizeof(int        ));
  PDM_g_num_t* delmt_ridge_edge_vtx     = (PDM_g_num_t *) malloc(  n_sum_vtx_ridge_edge_tot * sizeof(PDM_g_num_t));

  delmt_ridge_edge_vtx_idx[0] = 0;
  PDM_sections_decompose_edges(dmesh_nodal->ridge,
                               delmt_ridge_edge_vtx_idx,
                               delmt_ridge_edge_vtx,
                               delmt_ridge_edge,
                               NULL, NULL,
                               dparent_elmt_ridge_pos);

  /*
   *  Create empty dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     0,
                                     0,
                                     0,
                                     dmesh_nodal->vtx->n_vtx,
                                     dmesh_nodal->comm);

  /* Juste a view */
  // dm->_dvtx_coord = dmesh_nodal->vtx->_coords;
  PDM_dmesh_vtx_coord_set(dm, dmesh_nodal->vtx->_coords, PDM_OWNERSHIP_USER);

  assert(link->dmesh == NULL);
  link->dmesh = dm;

  /*
   * Tricks : We use the property of part_to_block to make the concatenation
   *      One part for surfacic
   *      Another for ridge
   */
  int n_part = 2;
  int* n_edge_elt_tot = malloc( n_part * sizeof(int));
  n_edge_elt_tot[0] = n_edge_elt_surf_tot;
  n_edge_elt_tot[1] = n_edge_elt_ridge_tot;

  PDM_g_num_t **delmt_edge         = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  int         **dparent_elmt_pos   = (int         **) malloc( n_part * sizeof(int         *));
  int         **delmt_edge_vtx_idx = (int         **) malloc( n_part * sizeof(int         *));
  PDM_g_num_t **delmt_edge_vtx     = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));

  delmt_edge        [0] = delmt_surf_edge        ;
  dparent_elmt_pos  [0] = dparent_elmt_surf_pos  ;
  delmt_edge_vtx_idx[0] = delmt_surf_edge_vtx_idx;
  delmt_edge_vtx    [0] = delmt_surf_edge_vtx    ;

  delmt_edge        [1] = delmt_ridge_edge        ;
  dparent_elmt_pos  [1] = dparent_elmt_ridge_pos  ;
  delmt_edge_vtx_idx[1] = delmt_ridge_edge_vtx_idx;
  delmt_edge_vtx    [1] = delmt_ridge_edge_vtx    ;


  if(dmesh_nodal->ridge->dparent_idx != NULL) {
    free(dmesh_nodal->ridge->dparent_idx);
  }
  if(dmesh_nodal->ridge->dparent_gnum != NULL) {
    free(dmesh_nodal->ridge->dparent_gnum);
  }
  if(dmesh_nodal->ridge->dparent_sign != NULL) {
    free(dmesh_nodal->ridge->dparent_sign);
  }
  if(dmesh_nodal->ridge->delmt_child_distrib != NULL) {
    free(dmesh_nodal->ridge->delmt_child_distrib);
  }

  /* Memory is deallocated inside */
  int n_section_child = dmesh_nodal->ridge->n_section;
  PDM_g_num_t n_g_child = dmesh_nodal->ridge->section_distribution[n_section_child];
  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_FACE] = PDM_TRUE;
  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
  PDM_generate_entitiy_connectivity(dmesh_nodal->comm,
                                    dmesh_nodal->n_vtx_abs,
                                    n_g_child,
                                    n_part,
                                    n_edge_elt_tot,
                                    delmt_edge,
                                    delmt_edge_vtx_idx,
                                    delmt_edge_vtx,
                                    dparent_elmt_pos,
                                    &dm->dn_edge,
                                    &dm->edge_distrib,
                                    &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                    &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                    &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                    &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                    &link->_dedge_parent_element_position,
                                    &dmesh_nodal->ridge->dparent_idx,
                                    &dmesh_nodal->ridge->dparent_gnum,
                                    &dmesh_nodal->ridge->dparent_sign,
                                    &dmesh_nodal->ridge->delmt_child_distrib,
                                    &link->distrib_missing_ridge,
                                    &link->dmissing_ridge_parent_g_num);
  dm->n_g_edge = dm->edge_distrib[dmesh_nodal->n_rank];

  int dn_ridge = (int) (dmesh_nodal->ridge->delmt_child_distrib[dmesh_nodal->i_rank+1] -
                        dmesh_nodal->ridge->delmt_child_distrib[dmesh_nodal->i_rank]);
  // PDM_log_trace_connectivity_long (dmesh_nodal->ridge->dparent_idx,
  //                                  dmesh_nodal->ridge->dparent_gnum,
  //                                  dn_ridge,
  //                                  "dmesh_nodal->ridge->dparent_gnum :");
  free(n_edge_elt_tot    );
  free(delmt_edge        );
  free(dparent_elmt_pos  );
  free(delmt_edge_vtx_idx);
  free(delmt_edge_vtx    );

  PDM_g_num_t *_dedge_vtx          = dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX];
  int         *_dedge_vtx_idx      = dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX];

  PDM_g_num_t* _dedge_face_tmp     = dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_FACE];
  int        * _dedge_face_idx_tmp = dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE];

   if (0) {
    PDM_log_trace_connectivity_long(_dedge_face_idx_tmp, _dedge_face_tmp, dm->dn_edge, "dedge_face_tmp : ");
  }

  // Post_treat
  if (link->distrib_missing_ridge[dmesh_nodal->n_rank] == 0 && post_treat_result == 1) {

    // PDM_log_trace_connectivity_long (_dedge_face_idx_tmp,
    //                                  _dedge_face_tmp,
    //                                  dm->dn_edge,
    //                                  "_dedge_face_tmp :");

    /* All children have a parent */
    PDM_g_num_t *dedge_face = (PDM_g_num_t *) malloc( 2 * dm->dn_edge * sizeof(PDM_g_num_t));
    int         *dflip_edge = (int         *) malloc(     dm->dn_edge * sizeof(int        ));
    for(int i_edge = 0; i_edge < dm->dn_edge; ++i_edge) {
      dflip_edge[i_edge] = 1;
      int beg = _dedge_face_idx_tmp[i_edge];
      int n_connect_face = _dedge_face_idx_tmp[i_edge+1] - beg;
      if(n_connect_face == 1) {
        dedge_face[2*i_edge  ] = PDM_ABS(_dedge_face_tmp[beg]); // Attention on peut être retourner !!!!
        dedge_face[2*i_edge+1] = 0;
      } else {
        assert(n_connect_face == 2);
        dedge_face[2*i_edge  ] = PDM_ABS(_dedge_face_tmp[beg  ]);
        dedge_face[2*i_edge+1] = PDM_ABS(_dedge_face_tmp[beg+1]);
      }

      // Flip if the edge is in the other sens
      if( PDM_SIGN(_dedge_face_tmp[beg]) == -1 ) {
        dflip_edge[i_edge] = -1;
        int beg_edge_vtx = _dedge_vtx_idx[i_edge  ];
        int end_edge_vtx = _dedge_vtx_idx[i_edge+1];

        PDM_g_num_t tmp_swap;
        tmp_swap = _dedge_vtx[beg_edge_vtx];
        _dedge_vtx[beg_edge_vtx  ] = _dedge_vtx[end_edge_vtx-1];
        _dedge_vtx[end_edge_vtx-1] = tmp_swap;
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_int(dflip_edge,  dm->dn_edge, "dflip_edge::");
    }

    /*
     * Actualize parent sign
     */
    //int dn_ridge = dmesh_nodal->ridge->delmt_child_distrib[dmesh_nodal->i_rank+1] - dmesh_nodal->ridge->delmt_child_distrib[dmesh_nodal->i_rank];
    PDM_block_to_part_t* btp = PDM_block_to_part_create(dm->edge_distrib,
                                                        (const PDM_g_num_t **) &dmesh_nodal->ridge->dparent_gnum,
                                                        &dn_ridge,
                                                        1,
                                                        dmesh_nodal->comm);

    int** tmp_pedge_flip;
    int stride_one = 1;
    PDM_block_to_part_exch(btp,
                            sizeof(int),
                            PDM_STRIDE_CST_INTERLACED,
                            &stride_one,
                            (void *  )   dflip_edge,
                            NULL,
                            (void ***)  &tmp_pedge_flip);
    int *pedge_flip = tmp_pedge_flip[0];

    for(int i = 0; i < dn_ridge; ++i) {
      dmesh_nodal->ridge->dparent_sign[i] *= pedge_flip[i];
    }

    free(pedge_flip);
    free(tmp_pedge_flip);
    free(dflip_edge);
    PDM_block_to_part_free(btp);


    free(_dedge_face_tmp);
    free(_dedge_face_idx_tmp);
    _dedge_face_tmp     = NULL;
    _dedge_face_idx_tmp = NULL;
    dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_FACE] = dedge_face;
    dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE] = NULL;

    if(0 == 1) {
      PDM_log_trace_array_long(dedge_face, 2 * dm->dn_edge, "edge_face::");
      PDM_log_trace_array_long(_dedge_vtx, _dedge_vtx_idx[dm->dn_edge], "_dedge_vtx::");
    }
  }

  // Not necessary because the algorithim keep the parent to find out boundary condition
  if (link->distrib_missing_ridge[dmesh_nodal->n_rank] == 0 && post_treat_result == 1) {
    assert(_dedge_face_tmp     == NULL);
    assert(_dedge_face_idx_tmp == NULL);
    PDM_setup_connectivity_idx(dm->dn_edge,
                               2,
                               dm->dconnectivity[PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                               &_dedge_face_idx_tmp,
                               &_dedge_face_tmp);
    for(int i = 0; i < dm->dn_edge; ++i) {
      int sgn = 1;
      for(int j = _dedge_face_idx_tmp[i]; j < _dedge_face_idx_tmp[i+1]; ++j){
        _dedge_face_tmp[j] = _dedge_face_tmp[j] * sgn;
        sgn *= -1;
      }
    }
  } else {
    assert(_dedge_face_tmp     != NULL);
    assert(_dedge_face_idx_tmp != NULL);
  }

  int is_signed = 1;
  assert(dm->face_distrib == NULL);
  dm->face_distrib = (PDM_g_num_t * ) malloc( (dmesh_nodal->n_rank + 1 ) * sizeof(PDM_g_num_t));
  dm->face_distrib[0] = -1;
  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
  PDM_dconnectivity_transpose(dmesh_nodal->comm,
                              dm->edge_distrib,
                              dm->face_distrib,
                              _dedge_face_idx_tmp,
                              _dedge_face_tmp,
                              is_signed,
                              &dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                              &dm->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE]);

  int dn_face = dm->face_distrib[dmesh_nodal->i_rank+1] - dm->face_distrib[dmesh_nodal->i_rank];
  dm->dn_face  = dn_face;
  dm->n_g_face = dm->face_distrib[dmesh_nodal->n_rank];

  if(0 == 1) {
    PDM_log_trace_array_long(dm->dconnectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE], dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][dn_face], "dface_edge :: ");
    PDM_log_trace_array_int(dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE], dm->dn_edge+1, "dedge_face_idx :: ");
    PDM_log_trace_array_long(dm->dconnectivity[PDM_CONNECTIVITY_TYPE_EDGE_FACE], dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE][dm->dn_edge], "dedge_face :: ");
    PDM_log_trace_connectivity_long(dm->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                    dm->dconnectivity[PDM_CONNECTIVITY_TYPE_EDGE_FACE],
                                    dm->dn_edge, "PDM_CONNECTIVITY_TYPE_EDGE_FACE :: ");
  }

  if (link->distrib_missing_ridge[dmesh_nodal->n_rank] == 0 && post_treat_result == 1) {
    free(_dedge_face_tmp);
    free(_dedge_face_idx_tmp);
  }
}

static
void
_generate_vtx_from_dmesh_nodal
(
  _pdm_link_dmesh_nodal_to_dmesh_t *link,
  int                               post_treat_result
)
{
  PDM_UNUSED(post_treat_result);

  PDM_dmesh_nodal_t* dmesh_nodal = link->dmesh_nodal;
  assert(dmesh_nodal->ridge    != NULL);

  int n_edge_elt_ridge_tot     = 0;
  int n_sum_vtx_ridge_edge_tot = 0;

  PDM_dmesh_nodal_elmts_decompose_edges_get_size(dmesh_nodal->ridge, &n_edge_elt_ridge_tot, &n_sum_vtx_ridge_edge_tot);

  /*
   * Decompose ridge
   */
  PDM_g_num_t* delmt_ridge_edge         = (PDM_g_num_t *) malloc(  n_edge_elt_ridge_tot     * sizeof(PDM_g_num_t));
  int*         dparent_elmt_ridge_pos   = (int         *) malloc(  n_edge_elt_ridge_tot     * sizeof(int        ));
  int*         delmt_ridge_edge_vtx_idx = (int         *) malloc( (n_edge_elt_ridge_tot +1) * sizeof(int        ));
  PDM_g_num_t* delmt_ridge_edge_vtx     = (PDM_g_num_t *) malloc(  n_sum_vtx_ridge_edge_tot * sizeof(PDM_g_num_t));

  delmt_ridge_edge_vtx_idx[0] = 0;
  PDM_sections_decompose_edges(dmesh_nodal->ridge,
                               delmt_ridge_edge_vtx_idx,
                               delmt_ridge_edge_vtx,
                               delmt_ridge_edge,
                               NULL, NULL,
                               dparent_elmt_ridge_pos);

  free(dparent_elmt_ridge_pos);
  free(delmt_ridge_edge);
  free(delmt_ridge_edge_vtx_idx);

  /*
   *  Create empty dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(link->owner,
                                     0,
                                     0,
                                     n_edge_elt_ridge_tot,
                                     dmesh_nodal->vtx->n_vtx,
                                     dmesh_nodal->comm);

  PDM_dmesh_vtx_coord_set(dm, dmesh_nodal->vtx->_coords, PDM_OWNERSHIP_USER);

  assert(link->dmesh == NULL);
  link->dmesh = dm;
  dm->dconnectivity_idx    [PDM_CONNECTIVITY_TYPE_EDGE_VTX] = NULL;
  dm->dconnectivity        [PDM_CONNECTIVITY_TYPE_EDGE_VTX] = delmt_ridge_edge_vtx;
  dm->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = PDM_TRUE;

}

static
void
_translate_element_group_to_entity
(
 PDM_MPI_Comm  comm,
 PDM_g_num_t  *entity_distrib,
 PDM_g_num_t  *dgroup_elmt,
 int          *dgroup_elmt_idx,
 int           n_group_elmt,
 int          *dchild_elt_parent_idx,
 PDM_g_num_t  *dchild_elt_parent,
 PDM_g_num_t **dentity_bound,
 int         **dentity_bound_idx
)
{
  /* dchild_elt_parent : Pour chaque enfant le numero parent de l'entité */
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  /*
   * Prepare exchange protocol
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity_distrib,
                               (const PDM_g_num_t **) &dgroup_elmt,
                                                      &dgroup_elmt_idx[n_group_elmt],
                                                      1,
                                                      comm);

  /*
   * Exchange
   */
  PDM_g_num_t** part_group_data;

  int dn_entity = (int) (entity_distrib[i_rank+1] - entity_distrib[i_rank]);
  int *block_stride = malloc (sizeof(int) * dn_entity);
  for (int i = 0; i < dn_entity; i++) {
    block_stride[i] = dchild_elt_parent_idx[i+1] - dchild_elt_parent_idx[i];
  }
  int **part_stride = NULL;

  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          block_stride,
             (void *  )   dchild_elt_parent,
                          &part_stride,
             (void ***)  &part_group_data);
  free (part_stride[0]);
  free (part_stride);
  free(block_stride);

  PDM_g_num_t* _part_group_data = part_group_data[0];

  *dentity_bound_idx = (int * ) malloc( (n_group_elmt+1) * sizeof(int) );
  int* _dentity_bound_idx = *dentity_bound_idx;

  for(int i_group = 0; i_group < n_group_elmt+1; ++i_group) {
    _dentity_bound_idx[i_group] = dgroup_elmt_idx[i_group];
  }

  *dentity_bound = _part_group_data;
  free(part_group_data);

  if(0 == 1) {
    PDM_log_trace_array_int (_dentity_bound_idx, n_group_elmt+1                  , "_dentity_bound_idx:: ");
    PDM_log_trace_array_long(_part_group_data  , _dentity_bound_idx[n_group_elmt], "_dentity_bound:: ");
  }
  PDM_block_to_part_free(btp);
}

static
void
_translate_element_group_to_faces
(
 _pdm_link_dmesh_nodal_to_dmesh_t* link
)
{
  PDM_dmesh_nodal_t *dmesh_nodal = link->dmesh_nodal;
  PDM_dmesh_t       *dm          = link->dmesh;

  PDM_g_num_t *dface_bound;
  int         *dface_bound_idx;

  if(dmesh_nodal->surfacic != NULL && dmesh_nodal->surfacic->n_group_elmt > 0) {
    _translate_element_group_to_entity(dmesh_nodal->comm,
                                       dmesh_nodal->surfacic->delmt_child_distrib,
                                       dmesh_nodal->surfacic->dgroup_elmt,
                                       dmesh_nodal->surfacic->dgroup_elmt_idx,
                                       dmesh_nodal->surfacic->n_group_elmt,
                                       dmesh_nodal->surfacic->dparent_idx,
                                       dmesh_nodal->surfacic->dparent_gnum,
                                       &dface_bound,
                                       &dface_bound_idx);

    dm->is_owner_bound[PDM_BOUND_TYPE_FACE] = PDM_TRUE;
    dm->dbound_idx    [PDM_BOUND_TYPE_FACE] = dface_bound_idx;
    dm->dbound        [PDM_BOUND_TYPE_FACE] = dface_bound;
    // dm->n_bnd                               = dmesh_nodal->surfacic->n_group_elmt; // TODO : TO REMOVE
    dm->n_group_bnd   [PDM_BOUND_TYPE_FACE] = dmesh_nodal->surfacic->n_group_elmt;
  }
  else {
    dm->is_owner_bound[PDM_BOUND_TYPE_FACE] = PDM_TRUE;
    dm->dbound_idx    [PDM_BOUND_TYPE_FACE] = PDM_array_zeros_int(1);
    dm->dbound        [PDM_BOUND_TYPE_FACE] = NULL;
    // dm->n_bnd                               = 0;
    dm->n_group_bnd   [PDM_BOUND_TYPE_FACE] = 0;
  }

  // Par recursion on peut avoir les group de vertex ou de edge


}

static
void
_translate_element_group_to_edges
(
 _pdm_link_dmesh_nodal_to_dmesh_t* link
)
{
  PDM_dmesh_nodal_t *dmesh_nodal = link->dmesh_nodal;
  PDM_dmesh_t       *dm          = link->dmesh;

  PDM_g_num_t *dedge_bound;
  int         *dedge_bound_idx;

  if(dmesh_nodal->ridge != NULL && dmesh_nodal->ridge->n_group_elmt > 0) {
    _translate_element_group_to_entity(dmesh_nodal->comm,
                                       dmesh_nodal->ridge->delmt_child_distrib,
                                       dmesh_nodal->ridge->dgroup_elmt,
                                       dmesh_nodal->ridge->dgroup_elmt_idx,
                                       dmesh_nodal->ridge->n_group_elmt,
                                       dmesh_nodal->ridge->dparent_idx,
                                       dmesh_nodal->ridge->dparent_gnum,
                                       &dedge_bound,
                                       &dedge_bound_idx);
    dm->is_owner_bound[PDM_BOUND_TYPE_EDGE] = PDM_TRUE;
    dm->dbound_idx    [PDM_BOUND_TYPE_EDGE] = dedge_bound_idx;
    dm->dbound        [PDM_BOUND_TYPE_EDGE] = dedge_bound;
    // dm->n_bnd                               = dmesh_nodal->ridge->n_group_elmt; // TODO : TO REMOVE
    dm->n_group_bnd   [PDM_BOUND_TYPE_EDGE] = dmesh_nodal->ridge->n_group_elmt;
  }
  else {
    dm->is_owner_bound[PDM_BOUND_TYPE_EDGE] = PDM_TRUE;
    dm->dbound_idx    [PDM_BOUND_TYPE_EDGE] = PDM_array_zeros_int(1);
    dm->dbound        [PDM_BOUND_TYPE_EDGE] = NULL;
    // dm->n_bnd                               = 0;
    dm->n_group_bnd   [PDM_BOUND_TYPE_EDGE] = 0;
  }
}


static
void
_translate_element_group_to_vtx
(
 _pdm_link_dmesh_nodal_to_dmesh_t* link
)
{
  PDM_dmesh_nodal_t *dmesh_nodal = link->dmesh_nodal;
  PDM_dmesh_t       *dm          = link->dmesh;

  if (dmesh_nodal->corner != NULL && dmesh_nodal->corner->n_group_elmt > 0) {
    // on a dgroup->elt et delt->vtx, on veut dgroup->vtx
    PDM_block_to_part_t *btp = PDM_block_to_part_create(dmesh_nodal->corner->sections_std[0]->distrib,
                                 (const PDM_g_num_t **) &dmesh_nodal->corner->dgroup_elmt,
                                                        &dmesh_nodal->corner->dgroup_elmt_idx[dmesh_nodal->corner->n_group_elmt],
                                                        1,
                                                        dmesh_nodal->comm);


    int one = 1;
    PDM_g_num_t **tmp_dcorner_elt = NULL;
    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &one,
                (void   *) dmesh_nodal->corner->sections_std[0]->_connec,
                           NULL,
                (void ***) &tmp_dcorner_elt);

    PDM_g_num_t *dcorner_elt = tmp_dcorner_elt[0];
    free(tmp_dcorner_elt);

    PDM_block_to_part_free(btp);

    int *dcorner_elt_idx = malloc(sizeof(int) * (dmesh_nodal->corner->n_group_elmt+1));
    memcpy(dcorner_elt_idx,
           dmesh_nodal->corner->dgroup_elmt_idx,
           sizeof(int) * (dmesh_nodal->corner->n_group_elmt+1));

    // _translate_element_group_to_entity(dmesh_nodal->comm,
    //                                    dmesh_nodal->corner->delmt_child_distrib,
    //                                    dmesh_nodal->corner->dgroup_elmt,
    //                                    dmesh_nodal->corner->dgroup_elmt_idx,
    //                                    dmesh_nodal->corner->n_group_elmt,
    //                                    dmesh_nodal->corner->dparent_idx,
    //                                    dmesh_nodal->corner->dparent_gnum,
    //                                    &dcorner_elt,
    //                                    &dcorner_elt_idx);
    dm->is_owner_bound[PDM_BOUND_TYPE_VTX] = PDM_TRUE;
    dm->dbound_idx    [PDM_BOUND_TYPE_VTX] = dcorner_elt_idx;
    dm->dbound        [PDM_BOUND_TYPE_VTX] = dcorner_elt;
  //   dm->n_bnd                              = dmesh_nodal->corner->n_group_elmt; // TODO : TO REMOVE
    dm->n_group_bnd   [PDM_BOUND_TYPE_VTX] = dmesh_nodal->corner->n_group_elmt;
  }
}


static
_pdm_link_dmesh_nodal_to_dmesh_t*
_link_dmesh_nodal_to_dmesh_init
(
 PDM_ownership_t owner
)
{
  _pdm_link_dmesh_nodal_to_dmesh_t* link = malloc( sizeof(_pdm_link_dmesh_nodal_to_dmesh_t));

  link->dmesh_nodal       = NULL;
  link->dmesh             = NULL;

  link->owner             = owner;
  link->results_is_getted = PDM_FALSE;

  link->dn_elmt           = 0;
  link->elmt_distrib      = NULL;

  link->_delmt_face       = NULL;
  link->_delmt_face_idx   = NULL;

  link->_dface_elmt                    = NULL;
  link->_dface_elmt_idx                = NULL;
  link->_dface_parent_element_position = NULL;

  link->_delmt_edge                    = NULL;
  link->_delmt_edge_idx                = NULL;

  link->_dedge_elmt       = NULL;
  link->_dedge_elmt_idx   = NULL;
  link->_dedge_parent_element_position = NULL;

  link->distrib_missing_ridge       = NULL;
  link->dmissing_ridge_parent_g_num = NULL;

  link->distrib_missing_surface       = NULL;
  link->dmissing_surface_parent_g_num = NULL;

  return link;
}


static
void
_link_dmesh_nodal_to_dmesh_free
(
 _pdm_link_dmesh_nodal_to_dmesh_t* link
)
{
  link->dmesh_nodal      = NULL; /* On a pas l'onwership de cette structure */

  if(( link->owner == PDM_OWNERSHIP_KEEP ) ||
     ( link->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !link->results_is_getted)){
    if(link->dmesh != NULL) {
      PDM_dmesh_free(link->dmesh);
    }
  }

  if(link->elmt_distrib != NULL) {
    free(link->elmt_distrib);
    link->elmt_distrib = NULL;
  }

  if(link->_delmt_face != NULL) {
    free(link->_delmt_face);
    link->_delmt_face = NULL;
  }

  if(link->_delmt_face_idx != NULL) {
    free(link->_delmt_face_idx);
    link->_delmt_face_idx = NULL;
  }

  if(link->_dface_elmt != NULL) {
    free(link->_dface_elmt);
    link->_dface_elmt = NULL;
  }

  if(link->_dface_elmt_idx != NULL) {
    free(link->_dface_elmt_idx);
    link->_dface_elmt_idx = NULL;
  }

  if(link->_dface_parent_element_position != NULL) {
    free(link->_dface_parent_element_position);
    link->_dface_parent_element_position = NULL;
  }

  if(link->_dedge_elmt != NULL) {
    free(link->_dedge_elmt);
    link->_dedge_elmt = NULL;
  }

  if(link->_dedge_elmt_idx != NULL) {
    free(link->_dedge_elmt_idx);
    link->_dedge_elmt_idx = NULL;
  }

  if(link->_dedge_parent_element_position != NULL) {
    free(link->_dedge_parent_element_position);
    link->_dedge_parent_element_position = NULL;
  }

  if (link->distrib_missing_ridge != NULL) {
    free (link->distrib_missing_ridge);
    link->distrib_missing_ridge = NULL;
  }

  if (link->dmissing_ridge_parent_g_num != NULL) {
    free (link->dmissing_ridge_parent_g_num);
    link->dmissing_ridge_parent_g_num = NULL;
  }

  if (link->distrib_missing_surface != NULL) {
    free (link->distrib_missing_surface);
    link->distrib_missing_surface = NULL;
  }

  if (link->dmissing_surface_parent_g_num != NULL) {
    free (link->dmissing_surface_parent_g_num);
    link->dmissing_surface_parent_g_num = NULL;
  }

  free(link);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_nodal_to_dmesh_t*
PDM_dmesh_nodal_to_dmesh_create
(
const int             n_mesh,
const PDM_MPI_Comm    comm,
const PDM_ownership_t owner
)
{
  PDM_dmesh_nodal_to_dmesh_t *dmesh_nodal_to_dm = (PDM_dmesh_nodal_to_dmesh_t *) malloc(sizeof(PDM_dmesh_nodal_to_dmesh_t));

  dmesh_nodal_to_dm->comm              = comm;
  dmesh_nodal_to_dm->owner             = owner;
  dmesh_nodal_to_dm->results_is_getted = PDM_FALSE;
  dmesh_nodal_to_dm->n_mesh            = n_mesh;
  dmesh_nodal_to_dm->link              = malloc( n_mesh * sizeof(_pdm_link_dmesh_nodal_to_dmesh_t*) );
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    dmesh_nodal_to_dm->link[i_mesh] = _link_dmesh_nodal_to_dmesh_init(owner);
  }
  dmesh_nodal_to_dm->post_treat_result = 1;

  return (PDM_dmesh_nodal_to_dmesh_t *) dmesh_nodal_to_dm;
}

void
PDM_dmesh_nodal_to_dmesh_set_post_treat_result
(
 PDM_dmesh_nodal_to_dmesh_t *dmesh_nodal_to_dm,
 int                         post_treat_result
)
{
  dmesh_nodal_to_dm->post_treat_result = post_treat_result;

}

/**
 * \brief  Add dmesh_nodal
 */
void
PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal
(
        PDM_dmesh_nodal_to_dmesh_t *dmesh_nodal_to_dm,
  const int                         i_mesh,
        PDM_dmesh_nodal_t          *dmesh_nodal
)
{
  dmesh_nodal_to_dm->link[i_mesh]->dmesh_nodal = dmesh_nodal;
}


void
PDM_dmesh_nodal_to_dmesh_get_dmesh
(
        PDM_dmesh_nodal_to_dmesh_t  *dmesh_nodal_to_dm,
  const int                          i_mesh,
        PDM_dmesh_t                **dm
)
{
  *dm = dmesh_nodal_to_dm->link[i_mesh]->dmesh;
  dmesh_nodal_to_dm->results_is_getted = PDM_TRUE;
}


/**
 * \brief  Free
 */
void
PDM_dmesh_nodal_to_dmesh_free
(
  PDM_dmesh_nodal_to_dmesh_t* dmesh_nodal_to_dm
)
{
  if(( dmesh_nodal_to_dm->owner == PDM_OWNERSHIP_KEEP ) ||
     ( dmesh_nodal_to_dm->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !dmesh_nodal_to_dm->results_is_getted)){

    for(int i_mesh = 0; i_mesh < dmesh_nodal_to_dm->n_mesh; ++i_mesh) {
      _link_dmesh_nodal_to_dmesh_free(dmesh_nodal_to_dm->link[i_mesh]);
    }
  } else {
    for(int i_mesh = 0; i_mesh < dmesh_nodal_to_dm->n_mesh; ++i_mesh) {
      _link_dmesh_nodal_to_dmesh_free(dmesh_nodal_to_dm->link[i_mesh]);
    }
  }

  free(dmesh_nodal_to_dm->link);

  free(dmesh_nodal_to_dm);

}

void
PDM_dmesh_nodal_to_dmesh_compute
(
        PDM_dmesh_nodal_to_dmesh_t                 *dmesh_nodal_to_dm,
  const PDM_dmesh_nodal_to_dmesh_transform_t        transform_kind,
  const PDM_dmesh_nodal_to_dmesh_translate_group_t  transform_group_kind
)
{
  double t1 = PDM_MPI_Wtime();
  for(int i_mesh = 0; i_mesh < dmesh_nodal_to_dm->n_mesh; ++i_mesh) {

    if(dmesh_nodal_to_dm->link[i_mesh]->dmesh_nodal->mesh_dimension == 2) {
      assert(transform_kind != PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE);
    }

    switch (transform_kind) {

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE:
        {
          _generate_faces_from_dmesh_nodal(dmesh_nodal_to_dm->link[i_mesh], dmesh_nodal_to_dm->post_treat_result);
        }
        break;

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE:
        {
          _generate_edges_from_dmesh_nodal(dmesh_nodal_to_dm->link[i_mesh], dmesh_nodal_to_dm->post_treat_result );
        }
        break;

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_VTX:
        {
          _generate_vtx_from_dmesh_nodal(dmesh_nodal_to_dm->link[i_mesh], dmesh_nodal_to_dm->post_treat_result );
        }
        break;
    }
  }

  // abort();
  // Boundary management
  for(int i_mesh = 0; i_mesh < dmesh_nodal_to_dm->n_mesh; ++i_mesh) {
    // if(dmesh_nodal_to_dm->link[i_mesh]->dmesh_nodal->dgroup_elmt != NULL) {
      switch (transform_group_kind) {
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_NONE:
          break;
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE:
          {
            assert(transform_kind == PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE);
            _translate_element_group_to_faces(dmesh_nodal_to_dm->link[i_mesh]);
            _translate_element_group_to_edges(dmesh_nodal_to_dm->link[i_mesh]);
            _translate_element_group_to_vtx  (dmesh_nodal_to_dm->link[i_mesh]);
          }
          break;
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE:
          {
            _translate_element_group_to_edges(dmesh_nodal_to_dm->link[i_mesh]);
            _translate_element_group_to_vtx  (dmesh_nodal_to_dm->link[i_mesh]);
            // PDM_error (__FILE__, __LINE__, 0, "PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE not implemented \n");
          }
          break;
        case PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_VTX:
          {
            // Do nothing all is translationg naturally
          }
          break;
      }
    // }
  }
  end_timer_and_print("PDM_dmesh_nodal_to_dmesh_compute", dmesh_nodal_to_dm->comm, t1);
}

void
PDM_generate_entitiy_connectivity_raw
(
PDM_MPI_Comm   comm,
PDM_g_num_t    n_vtx_abs,
int            n_entity_elt_tot,
PDM_g_num_t   *delmt_entity,
int           *delmt_entity_vtx_idx,
PDM_g_num_t   *delmt_entity_vtx,
int           *dn_entity,
PDM_g_num_t  **entity_distrib,
int          **dentity_vtx_idx,
PDM_g_num_t  **dentity_vtx,
int          **dentity_elmt_idx,
PDM_g_num_t  **dentity_elmt
)
{

  double t1 = PDM_MPI_Wtime();

  /*
   * We are now all information flatten - we only need to compute hash_keys for each entitys
   */
  PDM_g_num_t* ln_to_gn = (PDM_g_num_t*) malloc( n_entity_elt_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t key_mod = 4 * n_vtx_abs;

  // Option des clés
  _compute_keys(n_entity_elt_tot,
                delmt_entity_vtx_idx,
                delmt_entity_vtx,
                ln_to_gn,
                key_mod);

  if(0 == 1) {
    log_trace("n_entity_elt_tot = %i \n", n_entity_elt_tot);
    PDM_log_trace_array_long(ln_to_gn, n_entity_elt_tot , "ln_to_gn:: ");
    PDM_log_trace_array_int(delmt_entity_vtx_idx  , n_entity_elt_tot+1 , "delmt_entity_vtx_idx:: ");
    PDM_log_trace_array_long(delmt_entity_vtx, delmt_entity_vtx_idx[n_entity_elt_tot] , "delmt_entity_vtx:: ");
  }

  /*
   * Prepare exchange by computing stride
   */
  int* delmt_entity_vtx_n = (int        *) malloc( n_entity_elt_tot * sizeof(int        ));
  int* stride_one         = (int        *) malloc( n_entity_elt_tot * sizeof(int        ));
  for(int i_entity = 0; i_entity < n_entity_elt_tot; ++i_entity) {
    delmt_entity_vtx_n[i_entity] = delmt_entity_vtx_idx[i_entity+1] - delmt_entity_vtx_idx[i_entity];
    stride_one[i_entity]       = 1;
  }
  free(delmt_entity_vtx_idx);

  // PDM_log_trace_array_int(delmt_entity_vtx_n  , n_entity_elt_tot , "delmt_entity_vtx_n:: ");

  // for(int i = 0; i < n_entity_elt_tot; ++i) {
  //   log_trace("delmt_entity_vtx_n[%i] = %i \n", i, delmt_entity_vtx_n[i]);
  // }
  /*
   * Setup part_to_block to filter all keys
   */
  double* weight = (double *) malloc( n_entity_elt_tot * sizeof(double));
  for(int i = 0; i < n_entity_elt_tot; ++i) {
    weight[i] = 1.;
  }

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &ln_to_gn,
                                                      &weight,
                                                      &n_entity_elt_tot,
                                                      1,
                                                      comm);
  free(weight);
  /*
   * Exchange data
   */
  int*         blk_tot_entity_vtx_n = NULL;
  PDM_g_num_t* blk_tot_entity_vtx   = NULL;

  int blk_tot_entity_vtx_size = PDM_part_to_block_exch(         ptb,
                                                              sizeof(PDM_g_num_t),
                                                              PDM_STRIDE_VAR_INTERLACED,
                                                              -1,
                                                              &delmt_entity_vtx_n,
                                                    (void **) &delmt_entity_vtx,
                                                              &blk_tot_entity_vtx_n,
                                                    (void **) &blk_tot_entity_vtx);

  int* blk_n_entity_per_key = NULL;
  int* blk_entity_vtx_n     = NULL;
  int blk_entity_vtx_n_size = PDM_part_to_block_exch(         ptb,
                                                            sizeof(int),
                                                            PDM_STRIDE_VAR_INTERLACED,
                                                            -1,
                                                            &stride_one,
                                                  (void **) &delmt_entity_vtx_n,
                                                            &blk_n_entity_per_key,
                                                  (void **) &blk_entity_vtx_n);
  free(delmt_entity_vtx_n);

  int*         blk_elmt_entity_elmt_stri = NULL;
  PDM_g_num_t* blk_elmt_entity_elmt      = NULL;
  int blk_entity_elmt_size = PDM_part_to_block_exch(         ptb,
                                                           sizeof(PDM_g_num_t),
                                                           PDM_STRIDE_VAR_INTERLACED,
                                                           -1,
                                                           &stride_one,
                                                 (void **) &delmt_entity,
                                                           &blk_elmt_entity_elmt_stri,
                                                 (void **) &blk_elmt_entity_elmt);

  /*
   *  Get the size of the current process bloc
   */
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);

  if( 0 == 1 ) {
    PDM_log_trace_array_int(blk_tot_entity_vtx_n, blk_size             , "blk_tot_entity_vtx_n:: ");
    PDM_log_trace_array_long(blk_tot_entity_vtx , blk_tot_entity_vtx_size, "blk_tot_entity_vtx:: "  );

    PDM_log_trace_array_int(blk_n_entity_per_key, blk_size         , "blk_n_entity_per_key:: ");
    PDM_log_trace_array_int(blk_entity_vtx_n    , blk_entity_vtx_n_size, "blk_entity_vtx_n:: ");

    PDM_log_trace_array_long(blk_elmt_entity_elmt, blk_entity_elmt_size, "blk_elmt_entity_elmt:: ");
  }

  PDM_part_to_block_free(ptb);
  free(stride_one);
  free(blk_elmt_entity_elmt_stri); // Same as blk_n_entity_per_key
  free(blk_tot_entity_vtx_n);

  /*
   * Get the max number of vertex of entitys
   */
  int* blk_entity_vtx_idx  = (int        *) malloc( (blk_entity_vtx_n_size+1) * sizeof(int        ));
  int n_max_entity_per_key = 0;
  int n_tot_entity_per_key = 0;
  for(int i_entity = 0; i_entity < blk_size; ++i_entity) {
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, blk_n_entity_per_key[i_entity]);
    n_tot_entity_per_key += blk_n_entity_per_key[i_entity];
  }

  int n_max_vtx       = 0;
  blk_entity_vtx_idx[0] = 0;
  for(int i_entity = 0; i_entity < blk_entity_vtx_n_size; ++i_entity) {
    n_max_vtx          = PDM_MAX(n_max_vtx         , blk_entity_vtx_n    [i_entity]);
    blk_entity_vtx_idx[i_entity+1] = blk_entity_vtx_idx[i_entity] + blk_entity_vtx_n[i_entity];
  }

  // PDM_log_trace_array_int(blk_entity_vtx_idx, blk_entity_vtx_n_size, "blk_entity_vtx_idx:: ");
  /*
   * We need to identify each uniques entitys :
   *      - We have multiple packet to treat
   *      - The connectivity can be sorted in place
   *      - Multiple case can occur :
   *           - Alone entity normaly boundary
   *           - Multiple entitys
   *           - Same entity, we remove and replace by the first
   */
  PDM_g_num_t* loc_entity_vtx_1 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  PDM_g_num_t* loc_entity_vtx_2 = (PDM_g_num_t *) malloc(  n_max_vtx               * sizeof(PDM_g_num_t) );
  int*         already_treat    = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  int*         same_entity_idx  = (int         *) malloc( (n_max_entity_per_key+1) * sizeof(int        ) );
  int*         sens_entity      = (int         *) malloc(  n_max_entity_per_key    * sizeof(int        ) );
  PDM_g_num_t *tmp_parent       = (PDM_g_num_t *) malloc(n_max_entity_per_key      * sizeof(PDM_g_num_t) );
  int         *order            = (int         *) malloc(n_max_entity_per_key      * sizeof(int        ) );

  /*
   * Allocate Memory - entity_vtx - entity_elmt
   */
  *dentity_vtx      = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  blk_tot_entity_vtx_size  );
  *dentity_vtx_idx  = (int         *) malloc( sizeof(int        ) * (blk_entity_vtx_n_size+1 ));
  *dentity_elmt     = (PDM_g_num_t *) malloc( sizeof(PDM_g_num_t) *  n_tot_entity_per_key );
  *dentity_elmt_idx = (int         *) malloc( sizeof(int)         * (blk_entity_elmt_size+1)  );

  PDM_g_num_t *_dentity_vtx      = *dentity_vtx;
  int         *_dentity_vtx_idx  = *dentity_vtx_idx;
  PDM_g_num_t *_dentity_elmt     = *dentity_elmt;
  int         *_dentity_elmt_idx = *dentity_elmt_idx;

  // printf("blk_entity_elmt_size::%i\n", blk_entity_elmt_size);
  // printf("n_tot_entity_per_key::%i\n", n_tot_entity_per_key);


   /* Init global numbering */

  int i_abs_entity = 0;
  _dentity_vtx_idx[0] = 0;

  int idx = 0;
  int idx_entity_vtx = 0;
  _dentity_elmt_idx[0] = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {
    // printf(" --- Number of conflicting keys :: %i \n", blk_n_entity_per_key[i_key]);

    int n_conflict_entitys = blk_n_entity_per_key[i_key];

    for (int i = 0; i < n_conflict_entitys; i++) {
      tmp_parent[i] = blk_elmt_entity_elmt[idx+i];
      order[i] = i;
    }
    PDM_sort_long(tmp_parent, order, n_conflict_entitys);

    /* Reset */
    for(int j = 0; j < n_conflict_entitys; ++j) {
      already_treat[j] = -1;
    }

    /* Loop over all entitys in conflict */
    for(int idx_entity = 0; idx_entity < n_conflict_entitys; ++idx_entity) {
      int i_entity = order[idx_entity];
      // printf("Number of vtx on entitys %i :: %i with index [%i] \n", i_entity, blk_entity_vtx_n[idx+i_entity], idx+i_entity);

      int n_vtx_entity_1 = blk_entity_vtx_n  [idx+i_entity];
      int beg_1          = blk_entity_vtx_idx[idx+i_entity];
      int idx_next_same_entity = 0;
      sens_entity[idx_next_same_entity] = 1;
      same_entity_idx[idx_next_same_entity++] = i_entity;

      if(already_treat[i_entity] != 1) {

        PDM_g_num_t key_1 = 0;
        int idx_min_1 = -1;
        PDM_g_num_t min_1 = n_vtx_abs+1;
        for(int j = 0; j < n_vtx_entity_1; ++j) {
          loc_entity_vtx_1[j] = blk_tot_entity_vtx[beg_1+j];
          key_1 += loc_entity_vtx_1[j];
          if(loc_entity_vtx_1[j] < min_1) {
            min_1 = loc_entity_vtx_1[j];
            idx_min_1 = j;
          };
        }
        PDM_quick_sort_long(loc_entity_vtx_1, 0, n_vtx_entity_1-1);

        for(int idx_entity2 = 0; idx_entity2 < n_conflict_entitys; ++idx_entity2) {
          int i_entity_next = order[idx_entity2];

          if (i_entity_next == i_entity) {
            continue;
          }
          //printf("conflict : i_entity = %d, i_entity_next = %d...\n", i_entity, i_entity_next);
          if(already_treat[i_entity_next] == 1) {
            continue;
          }

          int n_vtx_entity_2 = blk_entity_vtx_n[idx+i_entity_next];

          if(n_vtx_entity_1 == n_vtx_entity_2 ) {

            int beg_2 = blk_entity_vtx_idx[idx+i_entity_next];
            PDM_g_num_t key_2 = 0;
            int idx_min_2 = -1;
            PDM_g_num_t min_2 = n_vtx_abs+1;
            for(int j = 0; j < n_vtx_entity_1; ++j) {
              loc_entity_vtx_2[j] = blk_tot_entity_vtx[beg_2+j];
              key_2 += loc_entity_vtx_2[j];
              if(loc_entity_vtx_2[j] < min_2) {
                min_2 = loc_entity_vtx_2[j];
                idx_min_2 = j;
              };
            }
            PDM_quick_sort_long(loc_entity_vtx_2, 0, n_vtx_entity_2-1);

            assert(key_1 % key_mod == key_2 % key_mod);

            int is_same_entity = 1;
            for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              if(loc_entity_vtx_1[i_vtx] != loc_entity_vtx_2[i_vtx]) {
                is_same_entity = -1;
                break;
              }
            }

            if(is_same_entity == 1 ){
              // printf("idx_min_1 = %i | idx_min_2 = %i \n", idx_min_1, idx_min_2);
              // printf(" test1 :: ");
              // for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              //   printf(" %i", (int)blk_tot_entity_vtx[beg_1+i_vtx]);
              // }
              // printf(" \n");
              // printf(" test2 :: ");
              // for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
              //   printf(" %i", (int)blk_tot_entity_vtx[beg_2+i_vtx]);
              // }
              // printf(" \n");

              if(n_vtx_entity_1 == 2) {
                PDM_g_num_t i1 = blk_tot_entity_vtx[beg_1  ];
                PDM_g_num_t i2 = blk_tot_entity_vtx[beg_1+1];

                PDM_g_num_t j1 = blk_tot_entity_vtx[beg_2];
                PDM_g_num_t j2 = blk_tot_entity_vtx[beg_2+1];
                if(i1 == j1) {
                  sens_entity[idx_next_same_entity] = 1;
                } else {
                  assert(i1 == j2);
                  assert(i2 == j1);
                  sens_entity[idx_next_same_entity] = -1;
                }

              } else {
                // Determine the sens
                PDM_g_num_t i1 = blk_tot_entity_vtx[beg_1 +  idx_min_1                    ];
                PDM_g_num_t i2 = blk_tot_entity_vtx[beg_1 + (idx_min_1+1) % n_vtx_entity_1];

                PDM_g_num_t j1 = blk_tot_entity_vtx[beg_2 +  idx_min_2                    ];
                PDM_g_num_t j2 = blk_tot_entity_vtx[beg_2 + (idx_min_2+1) % n_vtx_entity_1];
                // printf(" i1 = %i | i2 = %i | j1 = %i | j2 = %i\n", i1, i2, j1, j2);
                assert(i1 == j1); // Panic
                if(i2 != j2) {
                  sens_entity[idx_next_same_entity] = -1;
                  // printf(" idx3 = %i \n", (idx_min_1+n_vtx_entity_1) % n_vtx_entity_1);
                  PDM_g_num_t i3 = blk_tot_entity_vtx[beg_1 + (idx_min_1+n_vtx_entity_1-1) % n_vtx_entity_1];
                  // printf(" i1 = %i | i2 = %i | i3 = %i | j1 = %i | j2 = %i\n", i1, i2, i3, j1, j2);
                  assert(i3 == j2);
                } else {
                  sens_entity[idx_next_same_entity] = 1;
                  assert(i2 == j2);
                }
              }

              // Check if same sens
              same_entity_idx[idx_next_same_entity++] = i_entity_next;
            }
          } /* End if same number of vertex */
        } /* End loop next entity */

         // printf("[%i] - same_entity_idx::", idx_next_same_entity);
         // for(int ii = 0; ii < idx_next_same_entity; ++ii) {
         //   printf(" %i", same_entity_idx[ii]);
         // }
         // printf("\n");

        _dentity_vtx_idx[i_abs_entity+1] = _dentity_vtx_idx[i_abs_entity] + n_vtx_entity_1;
        for(int i_vtx = 0; i_vtx < n_vtx_entity_1; ++i_vtx) {
          _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
          // Ecriture à partir du minumun
          // _dentity_vtx[idx_entity_vtx++] = blk_tot_entity_vtx[beg_1+i_vtx];
        }

        _dentity_elmt_idx[i_abs_entity+1] = _dentity_elmt_idx[i_abs_entity];
        for(int i = 0; i < idx_next_same_entity; ++i) {
          int i_same_entity = same_entity_idx[i];
          int sign = sens_entity[i];
          int next_idx = _dentity_elmt_idx[i_abs_entity+1]++;
          // printf("[%i] _dentity_elmt[%i] = %i \n ", i_abs_entity,  next_idx, sign*blk_elmt_entity_elmt[idx+i_same_entity]);
          _dentity_elmt[next_idx] = sign*blk_elmt_entity_elmt[idx+i_same_entity];
          // _dentity_elmt[_dentity_elmt_idx[i_abs_entity+1]++] = sign*blk_elmt_entity_elmt[idx+i_same_entity];
          already_treat[i_same_entity] = 1;
        }
        i_abs_entity++;

      } /* End loop already treated */

    } /* End loop entity in conflict */

    idx += n_conflict_entitys;

  }

  // printf(" realloc dentity_elmt : %i --> %i \n", n_tot_entity_per_key, _dentity_elmt_idx[i_abs_entity]);
  *dentity_elmt = realloc(*dentity_elmt, sizeof(PDM_g_num_t) *  _dentity_elmt_idx[i_abs_entity] );
  _dentity_elmt = *dentity_elmt;

  /*
   * Free all unused structure
   */
  free(loc_entity_vtx_1);
  free(loc_entity_vtx_2);
  free(blk_entity_vtx_idx);
  free(already_treat);
  free(same_entity_idx);
  free(sens_entity);
  free(delmt_entity);
  free(delmt_entity_vtx);
  free(ln_to_gn);
  free(blk_tot_entity_vtx);
  free(blk_n_entity_per_key);
  free(blk_entity_vtx_n);
  free(blk_elmt_entity_elmt);
  free(tmp_parent);
  free(order);

  if( 0 == 1 ){
    printf("i_abs_entity::%i \n", i_abs_entity+1);
    PDM_log_trace_array_int(_dentity_vtx_idx, i_abs_entity+1                   , "_dentity_vtx_idx:: " );
    PDM_log_trace_array_long(_dentity_vtx   , _dentity_vtx_idx[i_abs_entity]   , "_dentity_vtx:: "     );
    PDM_log_trace_array_int(_dentity_elmt_idx, i_abs_entity+1                  , "_dentity_elmt_idx:: ");
    PDM_log_trace_array_long(_dentity_elmt    , _dentity_elmt_idx[i_abs_entity], "_dentity_elmt:: "    );
  }

  /*
   * Fill up structure
   */
  *dn_entity = i_abs_entity;
  int _dn_entity = *dn_entity;

  /*
   * Realloc
   */
  *dentity_vtx_idx = (int *        ) realloc(*dentity_vtx_idx, (_dn_entity + 1) * sizeof(int * ) );
  _dentity_vtx_idx = *dentity_vtx_idx;

  *dentity_vtx     = (PDM_g_num_t *) realloc(*dentity_vtx    , _dentity_vtx_idx[_dn_entity] * sizeof(PDM_g_num_t * ));
  _dentity_vtx     = *dentity_vtx;

  /*
   * Generate absolute numerotation of entitys
   */
  *entity_distrib = _make_absolute_entity_numbering(_dn_entity, comm);

  if (0) {
    end_timer_and_print("PDM_generate_entitiy_connectivity", comm, t1);
  }
}




void
PDM_dmesh_nodal_to_dmesh_get_missing
(
 PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
 const int                   i_mesh,
 PDM_geometry_kind_t         geom_kind,
 PDM_g_num_t               **distrib_missing,
 PDM_g_num_t               **dmissing_parent_g_num
 )
{
  if (geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
    *distrib_missing       = dmn_to_dm->link[i_mesh]->distrib_missing_ridge;
    *dmissing_parent_g_num = dmn_to_dm->link[i_mesh]->dmissing_ridge_parent_g_num;
  }
  else if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
    *distrib_missing       = dmn_to_dm->link[i_mesh]->distrib_missing_surface;
    *dmissing_parent_g_num = dmn_to_dm->link[i_mesh]->dmissing_surface_parent_g_num;
  }
  else {
    PDM_error (__FILE__, __LINE__, 0, "PDM_dmesh_nodal_to_dmesh_get_missing : invalid geom_kind %d\n", geom_kind);
  }
}
