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
#include "pdm_part_extension_algorithm.h"
#include "pdm_part_extension_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_domain_interface.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
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
void
_part_extension_3d
(
 PDM_part_extension_t *part_ext
)
{
  PDM_UNUSED(part_ext);
}

static
void
_part_extension_2d
(
 PDM_part_extension_t *part_ext
)
{
  PDM_UNUSED(part_ext);

}

static
void
_part_extension_1d
(
 PDM_part_extension_t *part_ext
)
{
  PDM_UNUSED(part_ext);

}

static
void
_concatenate_part_graph
(
 int           pn_entity1,
 int          *part1_to_part2_idx,
 int          *part1_to_part2_interface,
 int          *bentity1_entity2_n,
 PDM_g_num_t  *bentity1_entity2_gnum,
 int          *bentity1_entity2_triplet,
 int          *bentity1_entity2_interface_n,
 int          *bentity1_entity2_interface_tot_n,
 int          *bentity1_entity2_interface,
 int          *pnext_bentity1_entity2_n,
 PDM_g_num_t  *pnext_bentity1_entity2_gnum,
 int          *pnext_bentity1_entity2_triplet,
 int          *pnext_bentity1_entity2_interface_n,
 int          *pnext_bentity1_entity2_interface_tot_n,
 int          *pnext_bentity1_entity2_interface,
 int         **concat_bentity1_entity2_n_out,
 PDM_g_num_t **concat_bentity1_entity2_gnum_out,
 int         **concat_bentity1_entity2_triplet_out,
 int         **concat_bentity1_entity2_interface_n_out,
 int         **concat_bentity1_entity2_interface_tot_n_out,
 int         **concat_bentity1_entity2_interface_out
)
{
  PDM_UNUSED(part1_to_part2_interface);

  int         *concat_bentity1_entity2_n               = malloc(pn_entity1 * sizeof(int));;
  int         *concat_bentity1_entity2_interface_tot_n = malloc(pn_entity1 * sizeof(int));;

  int n_concat_entity1_entity2      = 0;
  int n_concat_entity1_entity2_itrf = 0;
  for(int i = 0; i < pn_entity1; ++i) {
    concat_bentity1_entity2_n              [i] = 0;
    concat_bentity1_entity2_interface_tot_n[i] = 0;

    /* From previous in layout of part2 */
    concat_bentity1_entity2_n              [i] += bentity1_entity2_n[i];
    concat_bentity1_entity2_interface_tot_n[i] += bentity1_entity2_interface_tot_n[i];

    /* From recv in layout part1_to_part2 */
    for(int idx = part1_to_part2_idx[i]/3; idx < part1_to_part2_idx[i+1]/3; ++idx) {
      concat_bentity1_entity2_n              [i] += pnext_bentity1_entity2_n              [idx];
      concat_bentity1_entity2_interface_tot_n[i] += pnext_bentity1_entity2_interface_tot_n[idx];
    }

    // Update global count
    n_concat_entity1_entity2      += concat_bentity1_entity2_n              [i];
    n_concat_entity1_entity2_itrf += concat_bentity1_entity2_interface_tot_n[i];
  }

  /* Allocate */
  PDM_g_num_t *concat_bentity1_entity2_gnum        = malloc(    n_concat_entity1_entity2      * sizeof(PDM_g_num_t));
  int         *concat_bentity1_entity2_triplet     = malloc(3 * n_concat_entity1_entity2      * sizeof(int        ));
  int         *concat_bentity1_entity2_interface_n = malloc(    n_concat_entity1_entity2      * sizeof(int        ));
  int         *concat_bentity1_entity2_interface   = malloc(    n_concat_entity1_entity2_itrf * sizeof(int        ));

  n_concat_entity1_entity2      = 0;
  n_concat_entity1_entity2_itrf = 0;
  int idx_readb      = 0;
  int idx_readb_itrf = 0;
  int idx_readp      = 0;
  int idx_readp_itrf = 0;
  for(int i = 0; i < pn_entity1; ++i) {

    /* From previous in layout of part2 */
    for(int k = 0; k < bentity1_entity2_n[i]; ++k) {
      concat_bentity1_entity2_gnum       [  n_concat_entity1_entity2  ] = bentity1_entity2_gnum       [  idx_readb  ];
      concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2  ] = bentity1_entity2_triplet    [3*idx_readb  ];
      concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2+1] = bentity1_entity2_triplet    [3*idx_readb+1];
      concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2+2] = bentity1_entity2_triplet    [3*idx_readb+2];
      concat_bentity1_entity2_interface_n[  n_concat_entity1_entity2  ] = bentity1_entity2_interface_n[  idx_readb  ];

      for(int p = 0; p <  bentity1_entity2_interface_n[  idx_readb  ]; ++p) {
        concat_bentity1_entity2_interface[n_concat_entity1_entity2_itrf++] = bentity1_entity2_interface[idx_readb_itrf++];
      }
      idx_readb++;
      n_concat_entity1_entity2++;
    }

    /* From recv in layout part1_to_part2 */
    for(int idx = part1_to_part2_idx[i]/3; idx < part1_to_part2_idx[i+1]/3; ++idx) {

      /* From previous in layout of part2 */
      for(int k = 0; k < pnext_bentity1_entity2_n[idx]; ++k) {
        concat_bentity1_entity2_gnum       [  n_concat_entity1_entity2  ] = pnext_bentity1_entity2_gnum       [  idx_readp  ];
        concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2  ] = pnext_bentity1_entity2_triplet    [3*idx_readp  ];
        concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2+1] = pnext_bentity1_entity2_triplet    [3*idx_readp+1];
        concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2+2] = pnext_bentity1_entity2_triplet    [3*idx_readp+2];
        concat_bentity1_entity2_interface_n[  n_concat_entity1_entity2  ] = pnext_bentity1_entity2_interface_n[  idx_readp  ];

        for(int p = 0; p <  pnext_bentity1_entity2_interface_n[  idx_readp  ]; ++p) {
          concat_bentity1_entity2_interface[n_concat_entity1_entity2_itrf++] = pnext_bentity1_entity2_interface[idx_readp_itrf++];
        }
        idx_readp++;
        n_concat_entity1_entity2++;
      }
    }
  }


  *concat_bentity1_entity2_n_out               = concat_bentity1_entity2_n;
  *concat_bentity1_entity2_gnum_out            = concat_bentity1_entity2_gnum;
  *concat_bentity1_entity2_triplet_out         = concat_bentity1_entity2_triplet;
  *concat_bentity1_entity2_interface_n_out     = concat_bentity1_entity2_interface_n;
  *concat_bentity1_entity2_interface_tot_n_out = concat_bentity1_entity2_interface_tot_n;
  *concat_bentity1_entity2_interface_out       = concat_bentity1_entity2_interface;

}

static
void
_dump_graph_info
(
 int           pn_entity1,
 int          *bentity1_entity2_n,
 PDM_g_num_t  *bentity1_entity2_gnum,
 int          *bentity1_entity2_triplet,
 int          *bentity1_entity2_interface_n,
 int          *bentity1_entity2_interface_tot_n,
 int          *bentity1_entity2_interface
)
{

  PDM_UNUSED(bentity1_entity2_interface);

  log_trace(" -------------------------- graph info -------------------------- \n");
  int idx_readb      = 0;
  int idx_readb_itrf = 0;
  for(int i = 0; i < pn_entity1; ++i) {
    log_trace("i_entity = %i - bentity1_entity2_n = %i - bentity1_entity2_interface_tot_n = %i\n", i, bentity1_entity2_n[i], bentity1_entity2_interface_tot_n[i]);
    /* From previous in layout of part2 */
    for(int k = 0; k < bentity1_entity2_n[i]; ++k) {
      // log_trace("\t bentity1_entity2_interface_n = %i \n", bentity1_entity2_interface_n[idx_readb]);
      log_trace("\t (gnum=%i, [%i/%i/%i] ", bentity1_entity2_gnum[idx_readb],
                                              bentity1_entity2_triplet[3*idx_readb],
                                              bentity1_entity2_triplet[3*idx_readb+1],
                                              bentity1_entity2_triplet[3*idx_readb+2]);
      if(bentity1_entity2_interface_n[idx_readb] > 0) {
        log_trace(" - Interface : ");
      }
      for(int p = 0; p < bentity1_entity2_interface_n[idx_readb]; ++p) {
        log_trace("%i ", bentity1_entity2_interface[idx_readb_itrf++]);
      }
      log_trace("\n");

      idx_readb++;
    }
  }


}



static
void
exchange_and_concat_part_graph
(
 PDM_part_to_part_t  *ptp,
 int                  n_part_tot,
 int                 *pn_entity1,
 int                **part1_to_part2_idx,
 int                **part1_to_part2_interface,
 int                **bentity1_entity2_n,
 PDM_g_num_t        **bentity1_entity2_gnum,
 int                **bentity1_entity2_triplet,
 int                **bentity1_entity2_interface_n,
 int                **bentity1_entity2_interface_tot_n,
 int                **bentity1_entity2_interface,
 int               ***pnext_bentity1_entity2_n_out,
 PDM_g_num_t       ***pnext_bentity1_entity2_gnum_out,
 int               ***pnext_bentity1_entity2_triplet_out,
 int               ***pnext_bentity1_entity2_interface_n_out,
 int               ***pnext_bentity1_entity2_interface_tot_n_out,
 int               ***pnext_bentity1_entity2_interface_out
)
{
  int exch_request = 0;
  int         **pnext_bentity1_entity2_n    = NULL;
  PDM_g_num_t **pnext_bentity1_entity2_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   bentity1_entity2_n,
                (const void **)  bentity1_entity2_gnum,
                                 &pnext_bentity1_entity2_n,
                    (void ***)   &pnext_bentity1_entity2_gnum,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pnext_bentity1_entity2_triplet = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 3 * sizeof(int),
                (const int **)   bentity1_entity2_n,
                (const void **)  bentity1_entity2_triplet,
                                 &pnext_bentity1_entity2_n,
                    (void ***)   &pnext_bentity1_entity2_triplet,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  /*
   * Exchange stride or interface
   */
  int **pnext_bentity1_entity2_interface_n = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(int),
                (const int **)   bentity1_entity2_n,
                (const void **)  bentity1_entity2_interface_n,
                                 &pnext_bentity1_entity2_n,
                    (void ***)   &pnext_bentity1_entity2_interface_n,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);


  int **pnext_bentity1_entity2_interface_tot_n = NULL;
  int **pnext_bentity1_entity2_interface       = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(int),
                (const int **)   bentity1_entity2_interface_tot_n,
                (const void **)  bentity1_entity2_interface,
                                 &pnext_bentity1_entity2_interface_tot_n,
                    (void ***)   &pnext_bentity1_entity2_interface,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  /*
   * Tout a été échanger avec l'interface courante, donc on rajoute à l'arrivé la provenace
   */
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    //

    int n_interface_tot = 0;
    for(int i = 0; i < pn_entity1[i_part]; ++i) {
      for(int idx = part1_to_part2_idx[i_part][i]/3; idx < part1_to_part2_idx[i_part][i+1]/3; ++idx) {
        n_interface_tot += pnext_bentity1_entity2_interface_tot_n[i_part][idx];
        if(part1_to_part2_interface[i_part][idx] != 0) {
          n_interface_tot += pnext_bentity1_entity2_n[i_part][idx];
        }
      }
    }

    int *pupdate_bentity1_entity2_interface = malloc(n_interface_tot * sizeof(int));

    n_interface_tot = 0;
    int idx_read      = 0;
    int idx_read_itrf = 0;
    for(int i = 0; i < pn_entity1[i_part]; ++i) {
      for(int idx = part1_to_part2_idx[i_part][i]/3; idx < part1_to_part2_idx[i_part][i+1]/3; ++idx) {


        for(int k = 0; k < pnext_bentity1_entity2_n[i_part][idx]; ++k) {


          for(int p = 0; p < pnext_bentity1_entity2_interface_n[i_part][idx_read]; ++p) {
            pupdate_bentity1_entity2_interface[n_interface_tot] = pnext_bentity1_entity2_interface[i_part][idx_read_itrf++];
            n_interface_tot++;
          }

          if(part1_to_part2_interface[i_part][idx] != 0) {
            pnext_bentity1_entity2_interface_tot_n[i_part][idx] += 1;
            pnext_bentity1_entity2_interface_n    [i_part][idx_read] += 1;
            pupdate_bentity1_entity2_interface[n_interface_tot++] = -part1_to_part2_interface[i_part][idx];
          }
          idx_read++;
        }

      }
    }

    free(pnext_bentity1_entity2_interface[i_part]);
    pnext_bentity1_entity2_interface[i_part] = pupdate_bentity1_entity2_interface;


  }



  *pnext_bentity1_entity2_n_out               = pnext_bentity1_entity2_n;
  *pnext_bentity1_entity2_gnum_out            = pnext_bentity1_entity2_gnum;
  *pnext_bentity1_entity2_triplet_out         = pnext_bentity1_entity2_triplet;
  *pnext_bentity1_entity2_interface_n_out     = pnext_bentity1_entity2_interface_n;
  *pnext_bentity1_entity2_interface_tot_n_out = pnext_bentity1_entity2_interface_tot_n;
  *pnext_bentity1_entity2_interface_out       = pnext_bentity1_entity2_interface;
}



static
void
_update_in_part2_layout
(
 int   pn_entity1,
 int  *part1_to_part2_idx,
 int  *pnext_bentity1_entity2_n,
 int  *pnext_bentity1_entity2_interface_tot_n,
 int **update_bentity1_entity2_n_out,
 int **update_bentity1_entity2_interface_tot_n_out
)
{
  int *update_bentity1_entity2_n               = malloc(pn_entity1 * sizeof(int));
  int *update_bentity1_entity2_interface_tot_n = malloc(pn_entity1 * sizeof(int));
  for(int i = 0; i < pn_entity1; ++i) {
    update_bentity1_entity2_n              [i] = 0;
    update_bentity1_entity2_interface_tot_n[i] = 0;
    for(int idx = part1_to_part2_idx[i]/3; idx < part1_to_part2_idx[i+1]/3; ++idx) {
      update_bentity1_entity2_n              [i] += pnext_bentity1_entity2_n              [idx];
      update_bentity1_entity2_interface_tot_n[i] += pnext_bentity1_entity2_interface_tot_n[idx];
    }
  }
  *update_bentity1_entity2_n_out               = update_bentity1_entity2_n;
  *update_bentity1_entity2_interface_tot_n_out = update_bentity1_entity2_interface_tot_n;
}



static
void
_recurse_and_filter
(
 PDM_part_to_part_t  *ptp,
 int                  n_part_tot,
 int                 *pn_entity1,
 int                **part1_to_part2_idx,
 int                **part1_to_part2_interface,
 int                **pentity1_entity2_idx,
 int                **pentity1_entity2,
 PDM_g_num_t        **pentity2_ln_to_gn,
 int                **pextract_entity2_n,
 PDM_g_num_t        **pextract_entity2_gnum,
 int                **pextract_entity2_triplet,
 PDM_MPI_Comm         comm
)
{

  PDM_UNUSED(pextract_entity2_n);
  PDM_UNUSED(pextract_entity2_gnum);
  PDM_UNUSED(pextract_entity2_triplet);


  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /*
   * We wan to build the complete graphe betwenn entity1 and entity2
   * But at each step we need to fix it with transformation
   * We need to build :
   *     - bentity1_entity2_n
   *     - bentity1_entity2_gnum
   */
  int         **bentity1_entity2_n               = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **bentity1_entity2_gnum            = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **bentity1_entity2_triplet         = malloc(n_part_tot * sizeof(int         *));
  int         **bentity1_entity2_interface_tot_n = malloc(n_part_tot * sizeof(int         *));
  int         **bentity1_entity2_interface       = malloc(n_part_tot * sizeof(int         *));
  int         **bentity1_entity2_interface_n     = malloc(n_part_tot * sizeof(int         *));

  /* Init */
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    bentity1_entity2_n              [i_part] = malloc(pn_entity1[i_part] * sizeof(int));
    bentity1_entity2_interface_tot_n[i_part] = malloc(pn_entity1[i_part] * sizeof(int));

    int         *_pentity1_entity2_idx             = pentity1_entity2_idx            [i_part];
    int         *_pentity1_entity2                 = pentity1_entity2                [i_part];
    PDM_g_num_t *_pentity2_ln_to_gn                = pentity2_ln_to_gn               [i_part];
    // int         *_part1_to_part2_idx               = part1_to_part2_idx              [i_part];
    int         *_bentity1_entity2_n               = bentity1_entity2_n              [i_part];
    int         *_bentity1_entity2_interface_tot_n = bentity1_entity2_interface_tot_n[i_part];

    int n_entity1_entity2      = 0;
    int n_entity1_entity2_itrf = 0;
    for(int i = 0; i < pn_entity1[i_part]; ++i) {
      _bentity1_entity2_n              [i] = 0;
      _bentity1_entity2_interface_tot_n[i] = 0;
      _bentity1_entity2_n[i] += _pentity1_entity2_idx[i+1] - _pentity1_entity2_idx[i];
      n_entity1_entity2 += _bentity1_entity2_n[i];
    }

    bentity1_entity2_gnum       [i_part] = malloc(    n_entity1_entity2      * sizeof(PDM_g_num_t));
    bentity1_entity2_triplet    [i_part] = malloc(3 * n_entity1_entity2      * sizeof(int        ));
    bentity1_entity2_interface_n[i_part] = malloc(    n_entity1_entity2      * sizeof(int        ));
    bentity1_entity2_interface  [i_part] = malloc(    n_entity1_entity2_itrf * sizeof(int        ));
    PDM_g_num_t *_bentity1_entity2_gnum        = bentity1_entity2_gnum       [i_part];
    int         *_bentity1_entity2_triplet     = bentity1_entity2_triplet    [i_part];
    int         *_bentity1_entity2_interface_n = bentity1_entity2_interface_n[i_part];
    int         *_bentity1_entity2_interface   = bentity1_entity2_interface  [i_part];

    n_entity1_entity2 = 0;
    for(int i = 0; i < pn_entity1[i_part]; ++i) {
      for(int j = _pentity1_entity2_idx[i]; j < _pentity1_entity2_idx[i+1]; ++j) {
        int i_entity2 = PDM_ABS(_pentity1_entity2[j])-1;
        _bentity1_entity2_gnum       [n_entity1_entity2] = _pentity2_ln_to_gn[i_entity2];
        _bentity1_entity2_interface_n[n_entity1_entity2] = 0;

        _bentity1_entity2_interface_tot_n[i] += _bentity1_entity2_interface_n[n_entity1_entity2];

        _bentity1_entity2_triplet    [3*n_entity1_entity2  ] = i_rank;
        _bentity1_entity2_triplet    [3*n_entity1_entity2+1] = i_part;
        _bentity1_entity2_triplet    [3*n_entity1_entity2+2] = i_entity2;

        n_entity1_entity2++;
      }
    }

    PDM_log_trace_array_int (_bentity1_entity2_n              , pn_entity1[i_part]    , "_bentity1_entity2_n           ::");
    PDM_log_trace_array_int (_bentity1_entity2_interface_tot_n, pn_entity1[i_part]    , "_bentity1_entity2_n           ::");
    PDM_log_trace_array_long(_bentity1_entity2_gnum           , n_entity1_entity2     , "_bentity1_entity2_gnum        ::");
    PDM_log_trace_array_int (_bentity1_entity2_triplet        , 3 * n_entity1_entity2 , "_bentity1_entity2_gnum        ::");
    PDM_log_trace_array_int (_bentity1_entity2_interface_n    , n_entity1_entity2     , "_bentity1_entity2_interface_n ::");
    PDM_log_trace_array_int (_bentity1_entity2_interface      , n_entity1_entity2_itrf, "_bentity1_entity2_interface   ::");
  }

  /*
   * Previous concat
   *
   */
  int         **prev_concat_bentity1_entity2_n               = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **prev_concat_bentity1_entity2_gnum            = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **prev_concat_bentity1_entity2_triplet         = malloc(n_part_tot * sizeof(int         *));
  int         **prev_concat_bentity1_entity2_interface_n     = malloc(n_part_tot * sizeof(int         *));
  int         **prev_concat_bentity1_entity2_interface_tot_n = malloc(n_part_tot * sizeof(int         *));
  int         **prev_concat_bentity1_entity2_interface       = malloc(n_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    prev_concat_bentity1_entity2_n              [i_part] = bentity1_entity2_n              [i_part];
    prev_concat_bentity1_entity2_gnum           [i_part] = bentity1_entity2_gnum           [i_part];
    prev_concat_bentity1_entity2_triplet        [i_part] = bentity1_entity2_triplet        [i_part];
    prev_concat_bentity1_entity2_interface_n    [i_part] = bentity1_entity2_interface_n    [i_part];
    prev_concat_bentity1_entity2_interface_tot_n[i_part] = bentity1_entity2_interface_tot_n[i_part];
    prev_concat_bentity1_entity2_interface      [i_part] = bentity1_entity2_interface      [i_part];
  }

  int         **concat_bentity1_entity2_n               = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **concat_bentity1_entity2_gnum            = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **concat_bentity1_entity2_triplet         = malloc(n_part_tot * sizeof(int         *));
  int         **concat_bentity1_entity2_interface_n     = malloc(n_part_tot * sizeof(int         *));
  int         **concat_bentity1_entity2_interface_tot_n = malloc(n_part_tot * sizeof(int         *));
  int         **concat_bentity1_entity2_interface       = malloc(n_part_tot * sizeof(int         *));

  for(int i_step = 0; i_step < 2; ++i_step)  {

    log_trace("------- STEP %i ------- \n", i_step);
    /*
     * Il faudrait refaire l'échange qu'on fait en gnum_come_from avec ce layout pour l'appeler en part2
     * On doit avoir le même resultats
     */
    int         **pnext_bentity1_entity2_n               = NULL;
    PDM_g_num_t **pnext_bentity1_entity2_gnum            = NULL;
    int         **pnext_bentity1_entity2_triplet         = NULL;
    int         **pnext_bentity1_entity2_interface_n     = NULL;
    int         **pnext_bentity1_entity2_interface_tot_n = NULL;
    int         **pnext_bentity1_entity2_interface       = NULL;

    for(int i_part = 0; i_part < n_part_tot; ++i_part) {
      PDM_log_trace_array_int (bentity1_entity2_n              [i_part], pn_entity1[i_part]    , "bentity1_entity2_n               ::");
      PDM_log_trace_array_int (bentity1_entity2_interface_tot_n[i_part], pn_entity1[i_part]    , "bentity1_entity2_interface_tot_n ::");
    }
    exchange_and_concat_part_graph(ptp,
                                   n_part_tot,
                                   pn_entity1,
                                   part1_to_part2_idx,
                                   part1_to_part2_interface,
                                   bentity1_entity2_n,
                                   bentity1_entity2_gnum,
                                   bentity1_entity2_triplet,
                                   bentity1_entity2_interface_n,
                                   bentity1_entity2_interface_tot_n,
                                   bentity1_entity2_interface,
                                   &pnext_bentity1_entity2_n,
                                   &pnext_bentity1_entity2_gnum,
                                   &pnext_bentity1_entity2_triplet,
                                   &pnext_bentity1_entity2_interface_n,
                                   &pnext_bentity1_entity2_interface_tot_n,
                                   &pnext_bentity1_entity2_interface);

    for(int i_part = 0; i_part < n_part_tot; ++i_part) {

      // Il faut rajouter l'interface number
      _concatenate_part_graph(pn_entity1[i_part],
                              part1_to_part2_idx                          [i_part],
                              part1_to_part2_interface                    [i_part],
                              prev_concat_bentity1_entity2_n              [i_part],
                              prev_concat_bentity1_entity2_gnum           [i_part],
                              prev_concat_bentity1_entity2_triplet        [i_part],
                              prev_concat_bentity1_entity2_interface_n    [i_part],
                              prev_concat_bentity1_entity2_interface_tot_n[i_part],
                              prev_concat_bentity1_entity2_interface      [i_part],
                              pnext_bentity1_entity2_n                    [i_part],
                              pnext_bentity1_entity2_gnum                 [i_part],
                              pnext_bentity1_entity2_triplet              [i_part],
                              pnext_bentity1_entity2_interface_n          [i_part],
                              pnext_bentity1_entity2_interface_tot_n      [i_part],
                              pnext_bentity1_entity2_interface            [i_part],
                              &concat_bentity1_entity2_n                  [i_part],
                              &concat_bentity1_entity2_gnum               [i_part],
                              &concat_bentity1_entity2_triplet            [i_part],
                              &concat_bentity1_entity2_interface_n        [i_part],
                              &concat_bentity1_entity2_interface_tot_n    [i_part],
                              &concat_bentity1_entity2_interface          [i_part]);

      int *update_bentity1_entity2_n               = NULL;
      int *update_bentity1_entity2_interface_tot_n = NULL;
      _update_in_part2_layout(pn_entity1                            [i_part],
                              part1_to_part2_idx                    [i_part],
                              pnext_bentity1_entity2_n              [i_part],
                              pnext_bentity1_entity2_interface_tot_n[i_part],
                              &update_bentity1_entity2_n,
                              &update_bentity1_entity2_interface_tot_n);

      free(pnext_bentity1_entity2_n              [i_part]);
      free(pnext_bentity1_entity2_interface_tot_n[i_part]);
      pnext_bentity1_entity2_n              [i_part] = update_bentity1_entity2_n;
      pnext_bentity1_entity2_interface_tot_n[i_part] = update_bentity1_entity2_interface_tot_n;

      /*
       * Preapre for next step
       */
      free(prev_concat_bentity1_entity2_n              [i_part]);
      free(prev_concat_bentity1_entity2_gnum           [i_part]);
      free(prev_concat_bentity1_entity2_triplet        [i_part]);
      free(prev_concat_bentity1_entity2_interface_n    [i_part]);
      free(prev_concat_bentity1_entity2_interface_tot_n[i_part]);
      free(prev_concat_bentity1_entity2_interface      [i_part]);

      prev_concat_bentity1_entity2_n              [i_part] = concat_bentity1_entity2_n              [i_part];
      prev_concat_bentity1_entity2_gnum           [i_part] = concat_bentity1_entity2_gnum           [i_part];
      prev_concat_bentity1_entity2_triplet        [i_part] = concat_bentity1_entity2_triplet        [i_part];
      prev_concat_bentity1_entity2_interface_n    [i_part] = concat_bentity1_entity2_interface_n    [i_part];
      prev_concat_bentity1_entity2_interface_tot_n[i_part] = concat_bentity1_entity2_interface_tot_n[i_part];
      prev_concat_bentity1_entity2_interface      [i_part] = concat_bentity1_entity2_interface      [i_part];

      concat_bentity1_entity2_n              [i_part] = NULL;
      concat_bentity1_entity2_gnum           [i_part] = NULL;
      concat_bentity1_entity2_triplet        [i_part] = NULL;
      concat_bentity1_entity2_interface_n    [i_part] = NULL;
      concat_bentity1_entity2_interface_tot_n[i_part] = NULL;
      concat_bentity1_entity2_interface      [i_part] = NULL;

      if(i_step > 0) {
        free(bentity1_entity2_n              [i_part]);
        free(bentity1_entity2_gnum           [i_part]);
        free(bentity1_entity2_triplet        [i_part]);
        free(bentity1_entity2_interface_n    [i_part]);
        free(bentity1_entity2_interface_tot_n[i_part]);
        free(bentity1_entity2_interface      [i_part]);
      }

      bentity1_entity2_n              [i_part] = pnext_bentity1_entity2_n              [i_part];
      bentity1_entity2_gnum           [i_part] = pnext_bentity1_entity2_gnum           [i_part];
      bentity1_entity2_triplet        [i_part] = pnext_bentity1_entity2_triplet        [i_part];
      bentity1_entity2_interface_n    [i_part] = pnext_bentity1_entity2_interface_n    [i_part];
      bentity1_entity2_interface_tot_n[i_part] = pnext_bentity1_entity2_interface_tot_n[i_part];
      bentity1_entity2_interface      [i_part] = pnext_bentity1_entity2_interface      [i_part];

      pnext_bentity1_entity2_n              [i_part] = NULL;
      pnext_bentity1_entity2_gnum           [i_part] = NULL;
      pnext_bentity1_entity2_triplet        [i_part] = NULL;
      pnext_bentity1_entity2_interface_n    [i_part] = NULL;
      pnext_bentity1_entity2_interface_tot_n[i_part] = NULL;
      pnext_bentity1_entity2_interface      [i_part] = NULL;


      _dump_graph_info(pn_entity1[i_part],
                       prev_concat_bentity1_entity2_n              [i_part],
                       prev_concat_bentity1_entity2_gnum           [i_part],
                       prev_concat_bentity1_entity2_triplet        [i_part],
                       prev_concat_bentity1_entity2_interface_n    [i_part],
                       prev_concat_bentity1_entity2_interface_tot_n[i_part],
                       prev_concat_bentity1_entity2_interface      [i_part]);

    }

    /*
     * Swap ptr
     */

    free(pnext_bentity1_entity2_n              );
    free(pnext_bentity1_entity2_gnum           );
    free(pnext_bentity1_entity2_triplet        );
    free(pnext_bentity1_entity2_interface_n    );
    free(pnext_bentity1_entity2_interface_tot_n);
    free(pnext_bentity1_entity2_interface      );

  }

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(prev_concat_bentity1_entity2_n              [i_part]);
    free(prev_concat_bentity1_entity2_gnum           [i_part]);
    free(prev_concat_bentity1_entity2_triplet        [i_part]);
    free(prev_concat_bentity1_entity2_interface_n    [i_part]);
    free(prev_concat_bentity1_entity2_interface_tot_n[i_part]);
    free(prev_concat_bentity1_entity2_interface      [i_part]);
  }
  free(prev_concat_bentity1_entity2_n              );
  free(prev_concat_bentity1_entity2_gnum           );
  free(prev_concat_bentity1_entity2_triplet        );
  free(prev_concat_bentity1_entity2_interface_n    );
  free(prev_concat_bentity1_entity2_interface_tot_n);
  free(prev_concat_bentity1_entity2_interface      );


  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {
  //   free(concat_bentity1_entity2_n              [i_part]);
  //   free(concat_bentity1_entity2_gnum           [i_part]);
  //   free(concat_bentity1_entity2_triplet        [i_part]);
  //   free(concat_bentity1_entity2_interface_n    [i_part]);
  //   free(concat_bentity1_entity2_interface_tot_n[i_part]);
  //   free(concat_bentity1_entity2_interface      [i_part]);
  // }
  free(concat_bentity1_entity2_n              );
  free(concat_bentity1_entity2_gnum           );
  free(concat_bentity1_entity2_triplet        );
  free(concat_bentity1_entity2_interface_n    );
  free(concat_bentity1_entity2_interface_tot_n);
  free(concat_bentity1_entity2_interface      );


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(bentity1_entity2_gnum           [i_part]);
    free(bentity1_entity2_n              [i_part]);
    free(bentity1_entity2_interface_n    [i_part]);
    free(bentity1_entity2_interface      [i_part]);
    free(bentity1_entity2_triplet        [i_part]);
    free(bentity1_entity2_interface_tot_n[i_part]);
  }
  free(bentity1_entity2_gnum           );
  free(bentity1_entity2_n              );
  free(bentity1_entity2_interface_n    );
  free(bentity1_entity2_interface      );
  free(bentity1_entity2_triplet        );
  free(bentity1_entity2_interface_tot_n);


  return;


  // _merge_graph(ptp,
  //              n_part_tot,
  //              pn_entity1,
  //              part1_to_part2_idx,
  //              part1_to_part2_interface,
  //              bentity1_entity2_n,
  //              bentity1_entity2_gnum,
  //              bentity1_entity2_triplet,
  //              bentity1_entity2_interface_n,
  //              bentity1_entity2_interface_tot_n,
  //              pnext_bentity1_entity2_n,
  //              pnext_bentity1_entity2_gnum,
  //              pnext_bentity1_entity2_triplet,
  //              pnext_bentity1_entity2_interface_n,
  //              pnext_bentity1_entity2_interface_tot_n,
  //              pnext_bentity1_entity2_interface,
  //              comm);


  // /*
  //  *  Free all recv
  //  */
  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {
  //   free(pnext_bentity1_entity2_n              [i_part]);
  //   free(pnext_bentity1_entity2_gnum           [i_part]);
  //   free(pnext_bentity1_entity2_triplet        [i_part]);
  //   free(pnext_bentity1_entity2_interface_n    [i_part]);
  //   free(pnext_bentity1_entity2_interface_tot_n[i_part]);
  //   free(pnext_bentity1_entity2_interface      [i_part]);
  // }
  // free(pnext_bentity1_entity2_n              );
  // free(pnext_bentity1_entity2_gnum           );
  // free(pnext_bentity1_entity2_triplet        );
  // free(pnext_bentity1_entity2_interface_n    );
  // free(pnext_bentity1_entity2_interface_tot_n);
  // free(pnext_bentity1_entity2_interface      );







  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {
  //   free(bentity1_entity2_gnum           [i_part]);
  //   free(bentity1_entity2_n              [i_part]);
  //   free(bentity1_entity2_interface_n    [i_part]);
  //   free(bentity1_entity2_interface      [i_part]);
  //   free(bentity1_entity2_triplet        [i_part]);
  //   free(bentity1_entity2_interface_tot_n[i_part]);
  // }
  // free(bentity1_entity2_gnum           );
  // free(bentity1_entity2_n              );
  // free(bentity1_entity2_interface_n    );
  // free(bentity1_entity2_interface      );
  // free(bentity1_entity2_triplet        );
  // free(bentity1_entity2_interface_tot_n);


  // /* Init step */
  // int         **pborder_entity1_to_pentity2_n = malloc(n_part_tot * sizeof(int         *));
  // int         **pborder_gnum_and_interface_n  = malloc(n_part_tot * sizeof(int         *));
  // PDM_g_num_t **pborder_gnum_and_interface    = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {

  //   // Compute index related of each entity1
  //   // int n_part1_to_part2 = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;

  //   pborder_entity1_to_pentity2_n[i_part] = malloc( pn_entity1[i_part] * sizeof(int        ));
  //   pborder_gnum_and_interface_n [i_part] = malloc( pn_entity1[i_part] * sizeof(int        ));

  //   int n_send_gnum_and_interface = 0;
  //   for(int i_entity = 0; i_entity < pn_entity1[i_part]; ++i_entity) {
  //     pborder_entity1_to_pentity2_n[i_part][i_entity] = 0;
  //     pborder_gnum_and_interface_n [i_part][i_entity] = 0;

  //     for(int idx = part1_to_part2_idx[i_part][i_entity]/3; idx < part1_to_part2_idx[i_part][i_entity+1]/3; ++idx) {
  //       pborder_entity1_to_pentity2_n[i_part][i_entity] +=     pextract_entity2_n[i_part][idx];
  //       pborder_gnum_and_interface_n [i_part][i_entity] += 3 * pextract_entity2_n[i_part][idx];
  //     }
  //     n_send_gnum_and_interface += pborder_gnum_and_interface_n [i_part][i_entity];
  //   }

  //   // pborder_entity1_to_pentity2_idx[i_part] = PDM_array_new_idx_from_sizes_int(pborder_entity1_to_pentity2_n[i_part], pn_entity1[i_part]);
  //   pborder_gnum_and_interface[i_part] = malloc(n_send_gnum_and_interface * sizeof(PDM_g_num_t));

  //   n_send_gnum_and_interface = 0;
  //   int idx_read = 0;
  //   for(int i_entity = 0; i_entity < pn_entity1[i_part]; ++i_entity) {
  //     for(int idx = part1_to_part2_idx[i_part][i_entity]/3; idx < part1_to_part2_idx[i_part][i_entity+1]/3; ++idx) {
  //       // At first step we have only simple interface
  //       for(int k = 0; k < pextract_entity2_n[i_part][idx]; ++k) {
  //         pborder_gnum_and_interface[i_part][n_send_gnum_and_interface++] = 2;
  //         pborder_gnum_and_interface[i_part][n_send_gnum_and_interface++] = pextract_entity2_gnum   [i_part][idx_read++];
  //         pborder_gnum_and_interface[i_part][n_send_gnum_and_interface++] = part1_to_part2_interface[i_part][idx  ];
  //       }
  //     }
  //   }

  //   if(1 == 0){
  //     PDM_log_trace_array_int(pborder_entity1_to_pentity2_n[i_part], pn_entity1[i_part], "pborder_entity1_to_pentity2_n[i_part] ::");
  //     PDM_log_trace_array_long(pborder_gnum_and_interface_n[i_part], pn_entity1[i_part], "pborder_gnum_and_interface_n[i_part] ::");
  //     PDM_log_trace_array_long(pborder_gnum_and_interface[i_part], n_send_gnum_and_interface, "pborder_gnum_and_interface[i_part] ::");
  //   }

  // }

  // int exch_request = 0;
  // int         **pnext_gnum_and_interface_n    = NULL;
  // PDM_g_num_t **pnext_gnum_and_interface_gnum = NULL;
  // PDM_part_to_part_reverse_iexch(ptp,
  //                                PDM_MPI_COMM_KIND_P2P,
  //                                PDM_STRIDE_VAR_INTERLACED,
  //                                PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
  //                                1,
  //                                sizeof(PDM_g_num_t),
  //               (const int **)   pborder_gnum_and_interface_n,
  //               (const void **)  pborder_gnum_and_interface,
  //                                &pnext_gnum_and_interface_n,
  //                   (void ***)   &pnext_gnum_and_interface_gnum,
  //                                &exch_request);
  // PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  // int         **pnext_extract_entity2_n    = NULL;
  // PDM_g_num_t **pnext_extract_entity2_gnum = NULL;
  // PDM_part_to_part_reverse_iexch(ptp,
  //                                PDM_MPI_COMM_KIND_P2P,
  //                                PDM_STRIDE_VAR_INTERLACED,
  //                                PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
  //                                1,
  //                                sizeof(PDM_g_num_t),
  //               (const int **)   pborder_entity1_to_pentity2_n,
  //               (const void **)  pextract_entity2_gnum,
  //                                &pnext_extract_entity2_n,
  //                   (void ***)   &pnext_extract_entity2_gnum,
  //                                &exch_request);
  // PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  // int **pnext_extract_entity2_triplet = NULL;
  // PDM_part_to_part_reverse_iexch(ptp,
  //                                PDM_MPI_COMM_KIND_P2P,
  //                                PDM_STRIDE_VAR_INTERLACED,
  //                                PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
  //                                1,
  //                                3 * sizeof(int),
  //               (const int  **)  pborder_entity1_to_pentity2_n,
  //               (const void **)  pextract_entity2_triplet,
  //                                &pnext_extract_entity2_n,
  //                   (void ***)   &pnext_extract_entity2_triplet,
  //                                &exch_request);
  // PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);


  // if(1 == 1) { // Usefull to know how many data is transfer
  //   for(int i_part = 0; i_part < n_part_tot; ++i_part) {

  //     int n_part1_to_part2 = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;
  //     PDM_log_trace_array_int(pnext_extract_entity2_n[i_part], n_part1_to_part2, "pnext_extract_entity2_n ::");
  //     int n_recv_tot  = 0;
  //     int n_recv_tot2 = 0;
  //     for(int i = 0; i < n_part1_to_part2; ++i ) {
  //       n_recv_tot  += pnext_extract_entity2_n   [i_part][i];
  //       n_recv_tot2 += pnext_gnum_and_interface_n[i_part][i];
  //     }

  //     PDM_log_trace_array_long(pnext_extract_entity2_gnum   [i_part], n_recv_tot , "pnext_extract_entity2_gnum    ::");
  //     PDM_log_trace_array_long(pnext_gnum_and_interface_gnum[i_part], n_recv_tot2, "pnext_gnum_and_interface_gnum ::");
  //   }
  // }


  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {

  //   int n_read = 0;
  //   for(int i_entity = 0; i_entity < pn_entity1[i_part]; ++i_entity) {

  //     log_trace("-------------------------------- i_entity = %i \n", i_entity);
  //     for(int idx = part1_to_part2_idx[i_part][i_entity]/3; idx < part1_to_part2_idx[i_part][i_entity+1]/3; ++idx) {

  //       log_trace("-------------------------------- part1_to_part2_interface[%i][%i] = %i \n", i_part, idx, part1_to_part2_interface[i_part][idx]);
  //       for(int k = 0; k < pnext_extract_entity2_n[i_part][idx]; ++k) {
  //         int n_data = pnext_gnum_and_interface_gnum[i_part][n_read];
  //         log_trace("\t n_data = %i / n_read = %i \n", n_data, n_read);


  //         PDM_g_num_t *lpnext_gnum_and_interface_gnum = &pnext_gnum_and_interface_gnum[i_part][n_read+1];

  //         PDM_log_trace_array_long(lpnext_gnum_and_interface_gnum, n_data, "");


  //         n_read += n_data+1;
  //       }
  //     }
  //   }
  // }



  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {
  //   int n_read = 0;
  //   for(int i_entity = 0; i_entity < pn_entity1[i_part]; ++i_entity) {

  //     for(int idx = part1_to_part2_idx[i_part][i_entity]/3; idx < part1_to_part2_idx[i_part][i_entity+1]/3; ++idx) {




  //     }
  //   }
  // }




  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {
  //   free(pborder_gnum_and_interface_n[i_part]);
  //   free(pborder_gnum_and_interface  [i_part]);
  // }
  // free(pborder_gnum_and_interface_n);
  // free(pborder_gnum_and_interface  );

  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {
  //   free(pnext_extract_entity2_n      [i_part]);
  //   free(pnext_extract_entity2_gnum   [i_part]);
  //   free(pnext_extract_entity2_triplet[i_part]);

  //   free(pnext_gnum_and_interface_n   [i_part]);
  //   free(pnext_gnum_and_interface_gnum[i_part]);
  //   free(pborder_entity1_to_pentity2_n[i_part]);

  // }
  // free(pnext_extract_entity2_n      );
  // free(pnext_extract_entity2_triplet);
  // free(pnext_gnum_and_interface_n   );
  // free(pnext_gnum_and_interface_gnum);
  // free(pnext_extract_entity2_gnum   );
  // free(pborder_entity1_to_pentity2_n);


}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*
 * Translate and post-treated link by interface with on entity to another
 *  Exemple of use :
 *     - entity1 to face
 *     - entity1 to edge
 *     - face to cell
 */

void
PDM_part_extension_interface_by_entity1_to_interface_by_entity2
(
  PDM_part_domain_interface_t  *pdi,
  PDM_bound_type_t              entity1_bound,
  int                           n_domain,
  PDM_g_num_t                  *shift_by_domain_entity2,
  int                          *n_part,
  int                         **pn_entity1_in,
  PDM_g_num_t                ***pentity1_ln_to_gn_in,
  int                        ***pentity1_hint_in,
  int                         **pn_entity2_in,
  PDM_g_num_t                ***pentity2_ln_to_gn_in,
  int                        ***pentity2_entity1_idx_in,
  int                        ***pentity2_entity1_in,
  int                         **pn_entity2_extented_out,
  PDM_g_num_t                ***pentity2_extented_ln_to_gn_out,
  int                        ***pentity2_extented_to_pentity2_idx_out,
  int                        ***pentity2_extented_to_pentity2_triplet_out,
  int                        ***pentity2_extented_to_pentity2_interface_out,
  PDM_MPI_Comm                  comm
)
{
  int n_part_tot = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    n_part_tot += n_part[i_domain];
  }

  int          *pn_entity1            = malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity1_ln_to_gn     = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity1_hint         = NULL;
  int          *pn_entity2            = malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_ln_to_gn     = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity2_entity1_idx  = malloc(n_part_tot * sizeof(int         *));
  int         **pentity2_entity1      = malloc(n_part_tot * sizeof(int         *));

  if(pentity1_hint_in != NULL) {
    pentity1_hint = malloc(n_part_tot * sizeof(int         *));
  }

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      pn_entity1          [ln_part_tot] = pn_entity1_in          [i_dom][i_part];
      pentity1_ln_to_gn   [ln_part_tot] = pentity1_ln_to_gn_in   [i_dom][i_part];
      pn_entity2          [ln_part_tot] = pn_entity2_in          [i_dom][i_part];
      pentity2_ln_to_gn   [ln_part_tot] = pentity2_ln_to_gn_in   [i_dom][i_part];
      pentity2_entity1_idx[ln_part_tot] = pentity2_entity1_idx_in[i_dom][i_part];
      pentity2_entity1    [ln_part_tot] = pentity2_entity1_in    [i_dom][i_part];

      if(pentity1_hint_in != NULL && pentity1_hint_in[i_dom] != NULL) {
        pentity1_hint[ln_part_tot] = pentity1_hint_in[i_dom][i_part];
      }

      ln_part_tot += 1;
    }
  }


  int         **pentity1_entity2_idx  = NULL;
  int         **pentity1_entity2      = NULL;
  PDM_part_connectivity_transpose(ln_part_tot,
                                  pn_entity2,
                                  pn_entity1,
                                  pentity2_entity1_idx,
                                  pentity2_entity1,
                                  &pentity1_entity2_idx,
                                  &pentity1_entity2);


  int          *pn_entity1_num             = NULL;
  int         **pentity1_num               = NULL;
  int         **pentity1_opp_location_idx  = NULL;
  int         **pentity1_opp_location      = NULL;
  int         **pentity1_opp_interface_idx = NULL;
  int         **pentity1_opp_interface     = NULL;
  int         **pentity1_opp_sens          = NULL;
  PDM_g_num_t **pentity1_opp_gnum          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         entity1_bound,
                                         pn_entity1,
                                         pentity1_ln_to_gn,
                                         &pn_entity1_num,
                                         &pentity1_num,
                                         &pentity1_opp_location_idx,
                                         &pentity1_opp_location,
                                         &pentity1_opp_interface_idx,
                                         &pentity1_opp_interface,
                                         &pentity1_opp_sens,
                                         &pentity1_opp_gnum);

  /*
   * Prepare creation of part_to_part
   */
  int         **part1_to_part2_idx              = malloc(n_part_tot * sizeof(int         *));
  int         **part1_to_part2_triplet_idx      = NULL; //malloc(n_part_tot * sizeof(int *));
  int         **part1_to_part2_triplet          = malloc(n_part_tot * sizeof(int         *));
  int         **part1_to_part2_interface        = malloc(n_part_tot * sizeof(int         *));

  int         **part1_to_part2_entity2_n           = malloc(n_part_tot * sizeof(int         *));
  int         **part1_to_part2_entity2_triplet     = malloc(n_part_tot * sizeof(int         *));

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_g_num_t *part_distribution = PDM_compute_entity_distribution(comm, n_part_tot);

  int **ppart_entity1_proc_idx = NULL;
  int **ppart_entity1_part_idx = NULL;
  int **ppart_entity1          = NULL; // (i_entity1, i_proc, i_part, i_entity1_opp)
  PDM_part_generate_entity_graph_comm(comm,
                                      part_distribution,
                                      NULL,
                                      n_part_tot,
                                      pn_entity1,
               (const PDM_g_num_t **) pentity1_ln_to_gn,
               (const int         **) pentity1_hint,
                                      &ppart_entity1_proc_idx,
                                      &ppart_entity1_part_idx,
                                      &ppart_entity1,
                                      NULL);

  free(part_distribution);

  int* n_part_g = malloc(n_domain * sizeof(int));
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);
  int n_g_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    n_g_part_tot += n_part_g[i_dom];
  }

  // Count
  int li_part = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    int *n_part_shift = (int *) malloc( (n_rank) * sizeof(int));
    PDM_MPI_Allgather(&li_part,
                      1,
                      PDM_MPI_INT,
                      n_part_shift,
                      1,
                      PDM_MPI_INT,
                      comm);

    if(0 == 1) {
      PDM_log_trace_array_int(n_part_shift, n_rank+1, "n_part_shift ::");
    }

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      part1_to_part2_idx[li_part] = malloc((pn_entity1[li_part] + 1) *sizeof(int));
      part1_to_part2_idx[li_part][0] = 0;

      int *part1_to_part2_n = PDM_array_zeros_int(pn_entity1[li_part]);

      int n_entity_bound = ppart_entity1_part_idx[i_part][n_g_part_tot];

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity = ppart_entity1[i_part][4*idx_entity]-1;
        part1_to_part2_n[i_entity] += 1;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_entity1_num[li_part]; ++idx_entity) {
        int i_entity = pentity1_num[li_part][idx_entity];
        int n_opp = pentity1_opp_location_idx[li_part][idx_entity+1] - pentity1_opp_location_idx[li_part][idx_entity];
        part1_to_part2_n[i_entity] += n_opp;
      }


      for(int i_entity = 0; i_entity < pn_entity1[li_part]; ++i_entity) {
        part1_to_part2_idx[li_part][i_entity+1] = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_n[i_entity] = 0;
      }

      int n_connect_tot = part1_to_part2_idx[li_part][pn_entity1[li_part]];
      part1_to_part2_triplet  [li_part] = malloc(n_connect_tot   * sizeof(int));
      part1_to_part2_interface[li_part] = malloc(n_connect_tot/3 * sizeof(int));

      // printf("n_connect_tot = %i \n", n_connect_tot);

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity     = ppart_entity1[i_part][4*idx_entity]-1;
        int i_proc_opp   = ppart_entity1[i_part][4*idx_entity+1];
        int i_part_opp   = ppart_entity1[i_part][4*idx_entity+2]-1;
        int i_entity_opp = ppart_entity1[i_part][4*idx_entity+3]-1;

        int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_triplet[li_part][idx_write  ] = i_proc_opp;
        part1_to_part2_triplet[li_part][idx_write+1] = i_part_opp + n_part_shift[i_proc_opp];
        part1_to_part2_triplet[li_part][idx_write+2] = i_entity_opp;

        idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
        part1_to_part2_interface[li_part][idx_write] = 0;

        // part1_to_part2_triplet[li_part][idx_entity+1] = part1_to_part2_triplet[li_part][idx_entity] + 3;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_entity1_num[li_part]; ++idx_entity) {
        int i_entity = pentity1_num[li_part][idx_entity];
        for(int idx_opp = pentity1_opp_location_idx[li_part][idx_entity  ];
                idx_opp < pentity1_opp_location_idx[li_part][idx_entity+1]; ++idx_opp) {

          int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
          part1_to_part2_triplet  [li_part][idx_write  ] = pentity1_opp_location [li_part][3*idx_opp  ];
          part1_to_part2_triplet  [li_part][idx_write+1] = pentity1_opp_location [li_part][3*idx_opp+1];
          part1_to_part2_triplet  [li_part][idx_write+2] = pentity1_opp_location [li_part][3*idx_opp+2];

          // Il faudra le faire en stride variable si periodicité composé
          idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
          part1_to_part2_interface[li_part][idx_write  ] = pentity1_opp_interface[li_part][  idx_opp  ];
        }
      }


      /*
       * Creation des buffers d'envoi des connectivités
       */
      int n_send = part1_to_part2_idx[li_part][pn_entity1[li_part]]/3;
      part1_to_part2_entity2_n    [li_part] = malloc(n_send * sizeof(int));
      int n_send_entity2     = 0;

      int *_pentity1_entity2_idx     = pentity1_entity2_idx    [li_part];
      int *_part1_to_part2_entity2_n = part1_to_part2_entity2_n[li_part];

      for(int i_entity = 0; i_entity < pn_entity1[li_part]; ++i_entity) {

        for(int idx = part1_to_part2_idx  [li_part][i_entity]/3; idx < part1_to_part2_idx  [li_part][i_entity+1]/3; ++idx) {
          _part1_to_part2_entity2_n    [idx] = 0;
          for(int idx_entity2 = _pentity1_entity2_idx[i_entity]; idx_entity2 < _pentity1_entity2_idx[i_entity+1]; ++idx_entity2 ) {
            // int i_entity2 = PDM_ABS(pentity1_entity2[li_part][idx_entity2])-1;
            _part1_to_part2_entity2_n    [idx] += 1;
          }
          n_send_entity2 += _part1_to_part2_entity2_n    [idx];
        }
      }


      part1_to_part2_entity2_triplet    [li_part] = malloc(3 * n_send_entity2     * sizeof(int        ));
      int         *_part1_to_part2_entity2_triplet   = part1_to_part2_entity2_triplet  [li_part];

      /* Fill loop */
      n_send_entity2     = 0;
      for(int i_entity = 0; i_entity < pn_entity1[li_part]; ++i_entity) {

        for(int idx = part1_to_part2_idx  [li_part][i_entity]/3; idx < part1_to_part2_idx  [li_part][i_entity+1]/3; ++idx) {

          /* Copy */
          for(int idx_entity2 = _pentity1_entity2_idx[i_entity]; idx_entity2 < _pentity1_entity2_idx[i_entity+1]; ++idx_entity2 ) {
            int i_entity2 = PDM_ABS(pentity1_entity2[li_part][idx_entity2])-1;

            _part1_to_part2_entity2_triplet[3*n_send_entity2  ] = i_rank;
            _part1_to_part2_entity2_triplet[3*n_send_entity2+1] = li_part;
            _part1_to_part2_entity2_triplet[3*n_send_entity2+2] = i_entity2;
            n_send_entity2++;

          }
        }
      }

      free(part1_to_part2_n);

      li_part += 1;
    }

    free(n_part_shift);
  }

  free(n_part_g);

  // TODO : Faire un part_to_block pour hook TOUTES les faces qui pointe sur le même sommet par exemple


  /*
   * Create part_to_part to exchange all data in opposit part
   */
  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pentity1_ln_to_gn,
                                                                      (const int          *) pn_entity1,
                                                                      n_part_tot,
                                                                      (const int          *) pn_entity1,
                                                                      n_part_tot,
                                                                      (const int         **) part1_to_part2_idx,
                                                                      (const int         **) part1_to_part2_triplet_idx,
                                                                      (const int         **) part1_to_part2_triplet,
                                                                      comm);


  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  /*
   * Prepare buffer
   */
  int         **gnum1_com_from_triplet_n    = malloc(n_part_tot * sizeof(int         *));
  int         **gnum1_com_from_triplet_send = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **gnum1_com_from_gnum_send    = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    int n_gnum1_come_from = gnum1_come_from_idx[i_part][ n_ref_lnum2[i_part]];
    int *_pentity1_entity2_idx     = pentity1_entity2_idx    [i_part];

    /* Count */
    int n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        n_send_part2 += _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
      }
    }

    /* Allocate */
    gnum1_com_from_triplet_n   [i_part] = malloc(    n_gnum1_come_from * sizeof(int        ));
    gnum1_com_from_triplet_send[i_part] = malloc(3 * n_send_part2      * sizeof(int        ));
    gnum1_com_from_gnum_send   [i_part] = malloc(    n_send_part2      * sizeof(PDM_g_num_t));

    int         *_gnum1_com_from_triplet_n    = gnum1_com_from_triplet_n   [i_part];
    int         *_gnum1_com_from_triplet_send = gnum1_com_from_triplet_send[i_part];
    PDM_g_num_t *_gnum1_com_from_gnum_send    = gnum1_com_from_gnum_send   [i_part];

    n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        _gnum1_com_from_triplet_n[k] = _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
        for(int idx_entity2 = _pentity1_entity2_idx[i_entity1]; idx_entity2 < _pentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {

          int i_entity2 = PDM_ABS(pentity1_entity2[i_part][idx_entity2])-1;
          _gnum1_com_from_gnum_send   [  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
          _gnum1_com_from_triplet_send[3*n_send_part2  ] = i_rank;
          _gnum1_com_from_triplet_send[3*n_send_part2+1] = i_part;
          _gnum1_com_from_triplet_send[3*n_send_part2+2] = i_entity2;
          n_send_part2++;
        }
      }
    }

    PDM_log_trace_array_long(_gnum1_com_from_gnum_send, n_send_part2, "_gnum1_com_from_gnum_send ::");


  }

  int exch_request = -1;
  int         **pextract_entity2_n    = NULL;
  PDM_g_num_t **pextract_entity2_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_gnum_send,
                                 &pextract_entity2_n,
                    (void ***)   &pextract_entity2_gnum,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pextract_entity2_triplet = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 3 * sizeof(int),
                (const int  **)  gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_triplet_send,
                                 &pextract_entity2_n,
                    (void ***)   &pextract_entity2_triplet,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(gnum1_com_from_triplet_n   [i_part]);
    free(gnum1_com_from_triplet_send[i_part]);
    free(gnum1_com_from_gnum_send   [i_part]);
  }


  if(1 == 1) { // Usefull to know how many data is transfer
    for(int i_part = 0; i_part < n_part_tot; ++i_part) {

      int n_triplet = part1_to_part2_idx[i_part][pn_entity1[i_part]];
      PDM_log_trace_array_int(part1_to_part2_idx       [i_part], pn_entity1[i_part]+1, "part1_to_part2_idx ::");
      PDM_log_trace_array_int(part1_to_part2_triplet   [i_part], n_triplet           , "part1_to_part2_triplet ::");
      PDM_log_trace_array_int(part1_to_part2_interface [i_part], n_triplet/3         , "part1_to_part2_interface ::");
      PDM_log_trace_graph_nuplet_int(part1_to_part2_idx[i_part],
                                     part1_to_part2_triplet[i_part], 1,
                                     pn_entity1[i_part], "part1_to_part2_triplet ::");


      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "ref_lnum2 :");

      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "extract_lnum2 :");
      PDM_log_trace_connectivity_long(gnum1_come_from_idx[i_part],
                                      gnum1_come_from    [i_part],
                                      n_ref_lnum2  [i_part], "gnum1_come_from ::");
    }
  }


  int **pextract_entity2_idx = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int n_part1_to_part2 = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;
    pextract_entity2_idx[i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity2_n[i_part], n_part1_to_part2);
  }


  _recurse_and_filter(ptp,
                      n_part_tot,
                      pn_entity1,
                      part1_to_part2_idx,
                      part1_to_part2_interface,
                      pentity1_entity2_idx,
                      pentity1_entity2,
                      pentity2_ln_to_gn,
                      pextract_entity2_n,
                      pextract_entity2_gnum,
                      pextract_entity2_triplet,
                      comm);

  free(gnum1_com_from_triplet_n   );
  free(gnum1_com_from_triplet_send);
  free(gnum1_com_from_gnum_send   );


  /*
   * Post-treatment
   *   - Remove duplicate (same gnum or same equivalent transformation )
   *   - Keep link between entity2 (extented and current part)
   *   - Create new global numbering for all new entities
   */
  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3,
                                                     n_part_tot,
                                                     PDM_TRUE,
                                                     1.e-6,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_entity2, 2);

  int          *pn_entity2_only_by_interface        = malloc(ln_part_tot * sizeof(int          ));
  int         **pentity2_interface                  = malloc(ln_part_tot * sizeof(int         *));
  PDM_g_num_t **pentity2_ln_to_gn_only_by_interface = malloc(ln_part_tot * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pextract_entity2_gnum    = pextract_entity2_gnum   [i_part];
    int         *_pextract_entity2_idx     = pextract_entity2_idx    [i_part];

    int n_part1_to_part2          = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot = _pextract_entity2_idx[n_part1_to_part2];
    pn_entity2_only_by_interface[i_part] = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(part1_to_part2_interface[i_part][i] != 0) {
        pn_entity2_only_by_interface[i_part] += pextract_entity2_n[i_part][i];
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_long(_pextract_entity2_gnum, n_part1_to_part2_recv_tot, "_pextract_entity2_gnum ::");
    }

    pentity2_ln_to_gn_only_by_interface[i_part] = malloc(2 * pn_entity2_only_by_interface[i_part] * sizeof(PDM_g_num_t));
    pentity2_interface                 [i_part] = malloc(    n_part1_to_part2_recv_tot            * sizeof(int        ));
    PDM_g_num_t *_pentity2_ln_to_gn_only_by_interface = pentity2_ln_to_gn_only_by_interface[i_part];

    pn_entity2_only_by_interface[i_part] = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(part1_to_part2_interface[i_part][i] != 0) {
        for(int j = _pextract_entity2_idx[i]; j < _pextract_entity2_idx[i+1]; ++j) {
          int idx_write = pn_entity2_only_by_interface[i_part]++;
          _pentity2_ln_to_gn_only_by_interface[2*idx_write  ] = _pextract_entity2_gnum[j];
          _pentity2_ln_to_gn_only_by_interface[2*idx_write+1] = part1_to_part2_interface[i_part][i];
        }
      }

      for(int j = _pextract_entity2_idx[i]; j < _pextract_entity2_idx[i+1]; ++j) {
        pentity2_interface[i_part][j] = part1_to_part2_interface[i_part][i];
      }
    }

    if(1 == 0) {
      PDM_log_trace_array_long(pentity2_ln_to_gn_only_by_interface[i_part], 2 * pn_entity2_only_by_interface[i_part], "pentity2_ln_to_gn_only_by_interface ::");
    }


    PDM_gnum_set_from_parents(gen_gnum_entity2,
                              i_part,
                              pn_entity2_only_by_interface[i_part],
                              pentity2_ln_to_gn_only_by_interface[i_part]);

  }

  PDM_gnum_compute(gen_gnum_entity2);

  /* Update */
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    PDM_g_num_t* extented_entity2_ln_to_gn = PDM_gnum_get(gen_gnum_entity2, i_part);

    PDM_g_num_t *_pextract_entity2_gnum    = pextract_entity2_gnum   [i_part];
    int         *_pextract_entity2_idx     = pextract_entity2_idx    [i_part];

    pn_entity2_only_by_interface[i_part] = 0;

    int n_part1_to_part2 = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;

    int idx_read = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(part1_to_part2_interface[i_part][i] != 0) {
        for(int j = _pextract_entity2_idx[i]; j < _pextract_entity2_idx[i+1]; ++j) {
          _pextract_entity2_gnum[j] = extented_entity2_ln_to_gn[idx_read++] + shift_by_domain_entity2[n_domain];
        }
      }
    }

    if(1 == 0) {
      PDM_log_trace_array_long(_pextract_entity2_gnum, _pextract_entity2_idx[n_part1_to_part2], "_pextract_entity2_gnum (Update) : ");
    }

  }


  PDM_gnum_free(gen_gnum_entity2);

  /*
   * Local post-treatment
   *   - Remove duplicate
   *   - Remove same entity but wih different path
   */
  int          *pn_entity2_extented                     = malloc(ln_part_tot * sizeof(int          ));
  int         **pentity2_extented_to_pentity2_idx       = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_to_pentity2_triplet   = malloc(ln_part_tot * sizeof(int         *));
  PDM_g_num_t **pentity2_extented_ln_to_gn              = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  PDM_g_num_t **extented_entity2_orig_gnum              = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity2_extented_to_pentity2_interface = malloc(ln_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pextract_entity2_gnum = pextract_entity2_gnum[i_part];
    int         *_pextract_entity2_idx  = pextract_entity2_idx [i_part];

    int n_part1_to_part2 = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot = _pextract_entity2_idx[n_part1_to_part2];

    // All gnum has be unified / shift w/r of interface, only sort along gnum
    int *order = malloc(n_part1_to_part2_recv_tot * sizeof(int));
    for(int i_entity2 = 0; i_entity2 < n_part1_to_part2_recv_tot; ++i_entity2) {
      order[i_entity2] = i_entity2;
    }

    int n_unique = PDM_inplace_unique_long_and_order(_pextract_entity2_gnum,
                                                     order,
                                                     0,
                                                     n_part1_to_part2_recv_tot-1);
    pn_entity2_extented[i_part] = n_unique;

    pentity2_extented_ln_to_gn[i_part] = malloc(n_unique * sizeof(PDM_g_num_t));
    for(int i = 0; i < n_unique; ++i) {
      pentity2_extented_ln_to_gn[i_part][i] = _pextract_entity2_gnum[i];
    }

    // Extract all data
    pentity2_extented_to_pentity2_idx      [i_part] = malloc( (     n_unique + 1) * sizeof(int        ));
    pentity2_extented_to_pentity2_triplet  [i_part] = malloc( ( 3 * n_unique    ) * sizeof(int        ));
    extented_entity2_orig_gnum             [i_part] = malloc( (     n_unique    ) * sizeof(PDM_g_num_t));
    pentity2_extented_to_pentity2_interface[i_part] = malloc( (     n_unique    ) * sizeof(int        ));

    /* Count */
    pentity2_extented_to_pentity2_idx    [i_part][0] = 0;
    for(int i_entity2 = 0; i_entity2 < n_unique; ++i_entity2) {
      int old_pos = order[i_entity2];
      pentity2_extented_to_pentity2_idx[i_part][i_entity2+1] = pentity2_extented_to_pentity2_idx[i_part][i_entity2] + 3;

      pentity2_extented_to_pentity2_triplet[i_part][3*i_entity2  ] = pextract_entity2_triplet[i_part][3*old_pos  ];
      pentity2_extented_to_pentity2_triplet[i_part][3*i_entity2+1] = pextract_entity2_triplet[i_part][3*old_pos+1];
      pentity2_extented_to_pentity2_triplet[i_part][3*i_entity2+2] = pextract_entity2_triplet[i_part][3*old_pos+2];

      // Save the link (gnum, interface)
      pentity2_extented_to_pentity2_interface[i_part][i_entity2] = pentity2_interface    [i_part][  old_pos  ];
      extented_entity2_orig_gnum             [i_part][i_entity2] = _pextract_entity2_gnum        [  i_entity2  ]; // Because it's sorted already
    }

    int n_triplet = pentity2_extented_to_pentity2_idx[i_part][n_unique];

    if(0 == 1) {
      PDM_log_trace_array_long(_pextract_entity2_gnum                         , n_unique   , "_pextract_entity2_gnum (UNIQUE)        ::");
      PDM_log_trace_array_int (pentity2_extented_to_pentity2_idx      [i_part], n_unique+1 , "pentity2_extented_to_pentity2_idx      ::");
      PDM_log_trace_array_long(extented_entity2_orig_gnum             [i_part], n_triplet/3, "extented_entity2_orig_gnum             ::");
      PDM_log_trace_array_int (pentity2_extented_to_pentity2_interface[i_part], n_triplet/3, "pentity2_extented_to_pentity2_interface::");
      PDM_log_trace_array_int (pentity2_extented_to_pentity2_triplet  [i_part], n_triplet  , "pentity2_extented_to_pentity2_triplet  ::");
    }

    free(order);

  }

  /*
   * A partir de toute ces infos on peut rajouter tout dans le PDM_part_domain_interface_set
   *   -> A généraliser pour des cellules
   *   -> Il manque les échanges pour les domaine
   *   -> Faire également un PDM_part_domain_interface_set_replace pour updater avec les nouvelles infos
   *   -> C'est cette info qu'on utilisera pour le rang 2 pour unifier les faces deja ramener
   */


  /*
   * Puis il faut également sortir toutes les infos pour construire le ptp pour faire la cascade
   *  part1 = part_extented
   *  part2 = part_current
   */
  *pn_entity2_extented_out                     = pn_entity2_extented;
  *pentity2_extented_ln_to_gn_out              = pentity2_extented_ln_to_gn;
  *pentity2_extented_to_pentity2_idx_out       = pentity2_extented_to_pentity2_idx;
  *pentity2_extented_to_pentity2_triplet_out   = pentity2_extented_to_pentity2_triplet;
  *pentity2_extented_to_pentity2_interface_out = pentity2_extented_to_pentity2_interface;



  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(extented_entity2_orig_gnum[i_part]);
  }
  free(extented_entity2_orig_gnum);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity2_ln_to_gn_only_by_interface[i_part]);
  }

  free(pn_entity2_only_by_interface);
  free(pentity2_ln_to_gn_only_by_interface);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_entity2_idx[i_part]);
    free(pentity2_interface  [i_part]);
  }
  free(pextract_entity2_idx);
  free(pentity2_interface);



  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_entity2_n        [i_part]);
    free(pextract_entity2_gnum     [i_part]);
    free(pextract_entity2_triplet  [i_part]);
    free(pentity1_entity2_idx      [i_part]);
    free(pentity1_entity2          [i_part]);
  }
  free(pextract_entity2_n        );
  free(pextract_entity2_gnum     );
  free(pextract_entity2_triplet  );
  free(pentity1_entity2_idx      );
  free(pentity1_entity2          );


  PDM_part_to_part_free(ptp);


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(ppart_entity1_proc_idx[i_part]);
    free(ppart_entity1_part_idx[i_part]);
    free(ppart_entity1         [i_part]);
  }
  free(ppart_entity1_proc_idx);
  free(ppart_entity1_part_idx);
  free(ppart_entity1         );


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(part1_to_part2_idx            [i_part] );
    free(part1_to_part2_triplet        [i_part]);
    free(part1_to_part2_interface      [i_part]);
    free(part1_to_part2_entity2_n      [i_part]);
    free(part1_to_part2_entity2_triplet[i_part]);

  }

  free(part1_to_part2_idx            );
  free(part1_to_part2_triplet        );
  free(part1_to_part2_interface      );
  free(part1_to_part2_entity2_n      );
  free(part1_to_part2_entity2_triplet);

  free(pn_entity1          );
  free(pentity1_ln_to_gn   );
  free(pn_entity2          );
  free(pentity2_ln_to_gn   );
  free(pentity2_entity1_idx);
  free(pentity2_entity1    );

  if(pentity1_hint_in != NULL) {
    free(pentity1_hint);
  }

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity1_num             [i_part]);
    free(pentity1_opp_location_idx[i_part]);
    free(pentity1_opp_location    [i_part]);
    free(pentity1_opp_interface   [i_part]);
    free(pentity1_opp_sens        [i_part]);
    free(pentity1_opp_gnum        [i_part]);
  }
  free(pn_entity1_num            );
  free(pentity1_num              );
  free(pentity1_opp_location_idx );
  free(pentity1_opp_location     );
  free(pentity1_opp_interface_idx);
  free(pentity1_opp_interface    );
  free(pentity1_opp_sens         );
  free(pentity1_opp_gnum         );

}




void
PDM_part_extension_pconnectivity_to_extented_pconnectivity
(
  PDM_part_domain_interface_t    *pdi,
  PDM_bound_type_t                entity2_bound,
  int                             n_domain,
  PDM_g_num_t                    *shift_by_domain_entity2,
  int                            *n_part,
  int                           **pn_entity1_in,
  PDM_g_num_t                  ***pentity1_ln_to_gn_in,
  int                           **pn_entity2_in,
  PDM_g_num_t                  ***pentity2_ln_to_gn_in,
  int                          ***pentity1_entity2_idx_in,
  int                          ***pentity1_entity2_in,
  int                            *pn_entity1_extented,
  PDM_g_num_t                   **pentity1_extented_ln_to_gn,
  int                           **pentity1_extented_to_pentity1_idx,
  int                           **pentity1_extented_to_pentity1_triplet,
  int                           **pentity1_extented_to_pentity1_interface,
  int                           **pn_entity2_extented_out,
  PDM_g_num_t                  ***pentity2_extented_ln_to_gn_out,
  int                          ***pextented_entity1_entity2_idx_out,
  int                          ***pextented_entity1_entity2_out,
  int                          ***pentity2_extented_to_pentity2_idx_out,
  int                          ***pentity2_extented_to_pentity2_triplet_out,
  int                          ***pentity2_extented_to_pentity2_interface_out,
  PDM_MPI_Comm                    comm
)
{
  log_trace("---  PDM_part_extension_pconnectivity_to_extented_pconnectivity --- \n");

  // Maybe supoose on entry already concat domain in part and juste give shift ?
  PDM_UNUSED(pdi);
  PDM_UNUSED(entity2_bound);
  PDM_UNUSED(n_domain);
  PDM_UNUSED(shift_by_domain_entity2);
  PDM_UNUSED(n_part);
  PDM_UNUSED(pn_entity1_in);
  PDM_UNUSED(pentity1_ln_to_gn_in);
  PDM_UNUSED(pn_entity2_in);
  PDM_UNUSED(pentity2_ln_to_gn_in);
  PDM_UNUSED(pentity1_entity2_idx_in);
  PDM_UNUSED(pentity1_entity2_in);
  PDM_UNUSED(pn_entity1_extented);
  PDM_UNUSED(pentity1_extented_ln_to_gn);
  PDM_UNUSED(pentity1_extented_to_pentity1_idx);
  PDM_UNUSED(pentity1_extented_to_pentity1_triplet);
  PDM_UNUSED(pentity1_extented_to_pentity1_interface);
  PDM_UNUSED(pn_entity2_extented_out);
  PDM_UNUSED(pentity2_extented_ln_to_gn_out);
  PDM_UNUSED(pextented_entity1_entity2_idx_out);
  PDM_UNUSED(pextented_entity1_entity2_out);
  PDM_UNUSED(pentity2_extented_to_pentity2_idx_out);
  PDM_UNUSED(pentity2_extented_to_pentity2_triplet_out);
  PDM_UNUSED(pentity2_extented_to_pentity2_interface_out);
  PDM_UNUSED(comm);

  PDM_domain_interface_t  *ditrf              = NULL;
  int                    **is_entity2_on_itrf = NULL;
  PDM_part_domain_interface_to_domain_interface(pdi,
                                                entity2_bound,
                                                n_part,
                                                pn_entity2_in,
                                                pentity2_ln_to_gn_in,
                                                &ditrf,
                                                &is_entity2_on_itrf);

  /*
   * Create flat view to have in flat frame all information for couple iterface / gnum
   */
  PDM_part_to_block_t **ptb_interface_entity2      = NULL;
  PDM_g_num_t         **interface_entity2_opp_gnum = NULL;
  PDM_domain_interface_make_flat_view(ditrf,
                                      entity2_bound,
                                      shift_by_domain_entity2,
                                      &ptb_interface_entity2,
                                      &interface_entity2_opp_gnum);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  int n_interface = PDM_part_domain_interface_n_interface_get(pdi);


  int n_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    n_part_tot += n_part[i_dom];
  }

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(is_entity2_on_itrf[i_part]);
  }
  free(is_entity2_on_itrf);

  PDM_domain_interface_free(ditrf);

  int          *pn_entity1           = (int          * ) malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity1_ln_to_gn    = (PDM_g_num_t ** ) malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity1_entity2_idx = (int         ** ) malloc(n_part_tot * sizeof(int         *));
  int         **pentity1_entity2     = (int         ** ) malloc(n_part_tot * sizeof(int         *));

  int          *pn_entity2        = (int          * ) malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_ln_to_gn = (PDM_g_num_t ** ) malloc(n_part_tot * sizeof(PDM_g_num_t *));

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      pn_entity1          [ln_part_tot+i_part] = pn_entity1_in          [i_dom][i_part];
      pentity1_ln_to_gn   [ln_part_tot+i_part] = pentity1_ln_to_gn_in   [i_dom][i_part];
      pentity1_entity2_idx[ln_part_tot+i_part] = pentity1_entity2_idx_in[i_dom][i_part];
      pentity1_entity2    [ln_part_tot+i_part] = pentity1_entity2_in    [i_dom][i_part];
      pn_entity2          [ln_part_tot+i_part] = pn_entity2_in          [i_dom][i_part];
      pentity2_ln_to_gn   [ln_part_tot+i_part] = pentity2_ln_to_gn_in   [i_dom][i_part];

    }
    ln_part_tot += n_part[i_dom];
  }

  /*
   * Preparation of array to manage interface at post-treatment
   */
  int          *pn_entity2_num             = NULL;
  int         **pentity2_num               = NULL;
  int         **pentity2_opp_location_idx  = NULL;
  int         **pentity2_opp_location      = NULL;
  int         **pentity2_opp_interface_idx = NULL;
  int         **pentity2_opp_interface     = NULL;
  int         **pentity2_opp_sens          = NULL;
  PDM_g_num_t **pentity2_opp_gnum          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         entity2_bound,
                                         pn_entity2,
                                         pentity2_ln_to_gn,
                                         &pn_entity2_num,
                                         &pentity2_num,
                                         &pentity2_opp_location_idx,
                                         &pentity2_opp_location,
                                         &pentity2_opp_interface_idx,
                                         &pentity2_opp_interface,
                                         &pentity2_opp_sens,
                                         &pentity2_opp_gnum);

  /*
   *  Il faut sortir le interface ln_to_gn + le shift puis garder un multiblock to part avec les doublé + les sens
   *
   */
  // Supprmier ce bout et utiliser le flat view APRès l'écahnge

  int          *pn_entity2_opp_gnum_and_itrf = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_opp_gnum_and_itrf   = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **entity2_opp_position         = malloc(ln_part_tot * sizeof(int         *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    // int *_pentity2_opp_location = pentity2_opp_location_idx[i_part];
    int *_pentity2_opp_location_idx = pentity2_opp_location_idx[i_part];
    int n_connect_tot = _pentity2_opp_location_idx[pn_entity2_num[i_part]];

    // On garde le lien courant (indirect sort of pentity2_num) + également keep le sens_opp
    pentity2_opp_gnum_and_itrf[i_part] = malloc(2 * n_connect_tot * sizeof(PDM_g_num_t));
    PDM_g_num_t* _pentity2_opp_gnum_and_itrf = pentity2_opp_gnum_and_itrf[i_part];

    // log_trace("pn_entity2_num[i_part] = %i \n", pn_entity2_num[i_part]);
    // PDM_log_trace_array_int(_pentity2_opp_location_idx, pn_entity2_num[i_part]+1, "_pentity2_opp_location_idx ::");

    int *tmp_opp_position = malloc(n_connect_tot * sizeof(int));
    for(int idx_entity = 0; idx_entity < pn_entity2_num[i_part]; ++idx_entity) {
      int i_entity = pentity2_num[i_part][idx_entity];
      for(int idx_opp = _pentity2_opp_location_idx[idx_entity]; idx_opp < _pentity2_opp_location_idx[idx_entity+1]; ++idx_opp) {

        _pentity2_opp_gnum_and_itrf[2*idx_opp  ] = pentity2_opp_gnum     [i_part][idx_opp];
        _pentity2_opp_gnum_and_itrf[2*idx_opp+1] = pentity2_opp_interface[i_part][idx_opp];

        tmp_opp_position[idx_opp] = i_entity;
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_long(pentity2_opp_gnum_and_itrf[i_part], 2 * n_connect_tot, "pentity2_opp_gnum_and_itrf ::");
      PDM_log_trace_array_int(tmp_opp_position, n_connect_tot, "tmp_opp_position ::");
    }

    entity2_opp_position        [i_part] = malloc(n_connect_tot * sizeof(int));

    int* order = malloc(n_connect_tot * sizeof(int));
    pn_entity2_opp_gnum_and_itrf[i_part] = PDM_order_inplace_unique_and_order_long(n_connect_tot,
                                                                                   2,
                                                                                   pentity2_opp_gnum_and_itrf[i_part],
                                                                                   order);


    if(0 == 1) {
      PDM_log_trace_array_int(order, n_connect_tot, "order = ");
      PDM_log_trace_array_long(pentity2_opp_gnum_and_itrf[i_part],
                               2 * pn_entity2_opp_gnum_and_itrf[i_part], "pentity2_opp_gnum_and_itrf ::");
      log_trace("pn_entity2_opp_gnum_and_itrf[i_part] = %i \n", pn_entity2_opp_gnum_and_itrf[i_part]);
      log_trace("n_connect_tot = %i \n", n_connect_tot);
    }

    for(int i = 0; i < pn_entity2_opp_gnum_and_itrf[i_part]; ++i) {
      entity2_opp_position[i_part][i] = tmp_opp_position[order[i]];
    }


    free(order);
    free(tmp_opp_position);

    if(0 == 1) {
      PDM_log_trace_array_int(entity2_opp_position[i_part],
                              pn_entity2_opp_gnum_and_itrf[i_part], "entity2_opp_position = ");
      PDM_log_trace_array_int(pentity2_num[i_part], pn_entity2_num[i_part], "pentity2_num ::");
    }
  }


  if(0 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      PDM_log_trace_array_int(pentity2_num[i_part], pn_entity2_num[i_part], "pentity2_num ::");
      PDM_log_trace_graph_nuplet_int(pentity2_opp_location_idx[i_part],
                                     pentity2_opp_location    [i_part],
                                     3,
                                     pn_entity2_num[i_part], "pentity2_opp_location ::");
      PDM_log_trace_array_int (pentity2_opp_interface[i_part], pentity2_opp_location_idx[i_part][pn_entity2_num[i_part]], "pentity2_opp_interface ::");
      PDM_log_trace_array_int (pentity2_opp_sens     [i_part], pentity2_opp_location_idx[i_part][pn_entity2_num[i_part]], "pentity2_opp_sens ::");
      PDM_log_trace_array_long(pentity2_opp_gnum     [i_part], pentity2_opp_location_idx[i_part][pn_entity2_num[i_part]], "pentity2_opp_gnum ::");
    }
  }

  /*
   * Creation du part_to_part entre entity2
   */
  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pentity1_extented_ln_to_gn,
                                                                      (const int          *) pn_entity1_extented,
                                                                      ln_part_tot,
                                                                      (const int          *) pn_entity1,
                                                                      ln_part_tot,
                                                                      (const int         **) pentity1_extented_to_pentity1_idx,
                                                                      (const int         **) NULL,
                                                                      (const int         **) pentity1_extented_to_pentity1_triplet,
                                                                      comm);

  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  /*
   * Prepare buffer
   */
  int         **gnum1_com_from_entity1_entity2_n       = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **gnum1_com_from_entity1_entity2         = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **gnum1_com_from_entity1_entity2_triplet = malloc(n_part_tot * sizeof(int         *));
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    int n_gnum1_come_from = gnum1_come_from_idx[i_part][ n_ref_lnum2[i_part]];
    int *_pentity1_entity2_idx = pentity1_entity2_idx[i_part];
    int *_pentity1_entity2     = pentity1_entity2    [i_part];

    /* Count */
    int n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        n_send_part2 += _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
      }
    }

    /* Allocate */
    gnum1_com_from_entity1_entity2_n      [i_part] = malloc(     n_gnum1_come_from * sizeof(int        ));
    gnum1_com_from_entity1_entity2        [i_part] = malloc(     n_send_part2      * sizeof(PDM_g_num_t));
    gnum1_com_from_entity1_entity2_triplet[i_part] = malloc( 3 * n_send_part2      * sizeof(int        ));

    int         *_gnum1_com_from_triplet_n               = gnum1_com_from_entity1_entity2_n      [i_part];
    PDM_g_num_t *_gnum1_com_from_entity1_entity2         = gnum1_com_from_entity1_entity2        [i_part];
    int         *_gnum1_com_from_entity1_entity2_triplet = gnum1_com_from_entity1_entity2_triplet[i_part];

    n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        _gnum1_com_from_triplet_n[k] = _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
        for(int idx_entity2 = _pentity1_entity2_idx[i_entity1]; idx_entity2 < _pentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {
          int i_entity2 = PDM_ABS(_pentity1_entity2[idx_entity2])-1;
          _gnum1_com_from_entity1_entity2        [  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
          _gnum1_com_from_entity1_entity2_triplet[3*n_send_part2  ] = i_rank;
          _gnum1_com_from_entity1_entity2_triplet[3*n_send_part2+1] = i_part;
          _gnum1_com_from_entity1_entity2_triplet[3*n_send_part2+2] = i_entity2;
          n_send_part2++;
        }
      }
    }
  }


  int exch_request = -1;
  int         **pextract_entity1_entity2_n    = NULL;
  PDM_g_num_t **pextract_entity1_entity2_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   gnum1_com_from_entity1_entity2_n,
                (const void **)  gnum1_com_from_entity1_entity2,
                                 &pextract_entity1_entity2_n,
                    (void ***)   &pextract_entity1_entity2_gnum,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pextract_entity1_entity2_triplet = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 3 * sizeof(int),
                (const int **)   gnum1_com_from_entity1_entity2_n,
                (const void **)  gnum1_com_from_entity1_entity2_triplet,
                                 &pextract_entity1_entity2_n,
                    (void ***)   &pextract_entity1_entity2_triplet,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(gnum1_com_from_entity1_entity2_n      [i_part]);
    free(gnum1_com_from_entity1_entity2        [i_part]);
    free(gnum1_com_from_entity1_entity2_triplet[i_part]);
  }
  free(gnum1_com_from_entity1_entity2_n      );
  free(gnum1_com_from_entity1_entity2        );
  free(gnum1_com_from_entity1_entity2_triplet);

  PDM_part_to_part_free(ptp);

  int **pextract_entity1_entity2_idx = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
    pextract_entity1_entity2_idx[i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity1_entity2_n[i_part], n_part1_to_part2);
  }


  /*
   * Get for all receive connectivity information about sens and opp gnum
   * CAUTION we need to manage SENS and SIGN
   * We concatenate for all partition
   */
  PDM_g_num_t ***query_itrf_gnum    = malloc(n_interface * sizeof(PDM_g_num_t **));
  int          **query_itrf_n       = malloc(n_interface * sizeof(int          *));
  int         ***recv_itrf_opp_n    = malloc(n_interface * sizeof(int         **));
  PDM_g_num_t ***recv_itrf_opp_gnum = malloc(n_interface * sizeof(PDM_g_num_t **));

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    query_itrf_gnum   [i_interface] = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
    recv_itrf_opp_n   [i_interface] = malloc(ln_part_tot * sizeof(int         *));
    recv_itrf_opp_gnum[i_interface] = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
    query_itrf_n      [i_interface] = malloc(ln_part_tot * sizeof(PDM_g_num_t  ));

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      query_itrf_n      [i_interface][i_part] = 0;
    }
  }


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int         *_pextract_entity1_interface    = pentity1_extented_to_pentity1_interface[i_part];
    int         *_pextract_entity1_entity2_idx  = pextract_entity1_entity2_idx    [i_part];
    PDM_g_num_t *_pextract_entity1_entity2_gnum = pextract_entity1_entity2_gnum   [i_part];

    /* Count */
    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
      if(_pextract_entity1_interface[i_ref] != 0) {
        int i_itrf = PDM_ABS(_pextract_entity1_interface[i_ref])-1;
        query_itrf_n[i_itrf][i_part] += _pextract_entity1_entity2_idx[i_ref+1]-_pextract_entity1_entity2_idx[i_ref];
      }
    }

    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      query_itrf_gnum   [i_interface][i_part] = malloc(query_itrf_n[i_interface][i_part] * sizeof(PDM_g_num_t));
      recv_itrf_opp_n   [i_interface][i_part] = malloc(query_itrf_n[i_interface][i_part] * sizeof(int        ));
      recv_itrf_opp_gnum[i_interface][i_part] = malloc(query_itrf_n[i_interface][i_part] * sizeof(PDM_g_num_t));
      query_itrf_n[i_interface][i_part] = 0;
    }

    /* Fill */
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
      if(_pextract_entity1_interface[i_ref] != 0) {
        int i_itrf = PDM_ABS(_pextract_entity1_interface[i_ref])-1;
        /* Blinder d'assert */
        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {
          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          int idx_write = query_itrf_n[i_itrf][i_part]++;
          query_itrf_gnum[i_itrf][i_part][idx_write] = entity2_g_num;
        }
      }
    }
  }


  if(1 == 1) {
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
        PDM_log_trace_array_long(query_itrf_gnum[i_interface][i_part], query_itrf_n[i_interface][i_part], "query_itrf_gnum ::");
      }
    }
  }

  PDM_g_num_t ***recv_interface_entity2_opp_gnum   = malloc(n_interface * sizeof(PDM_g_num_t **));
  int          **n_recv_interface_entity2_opp_gnum = malloc(n_interface * sizeof(int          *));
  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

    n_recv_interface_entity2_opp_gnum[i_interface] = malloc(ln_part_tot * sizeof(int));

    int n_gnum            = PDM_part_to_block_n_elt_block_get(ptb_interface_entity2[i_interface]);
    PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get (ptb_interface_entity2[i_interface]);

    PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block(blk_gnum,
                                                                          n_gnum,
                                                  (const PDM_g_num_t **)  query_itrf_gnum[i_interface],
                                                                          query_itrf_n   [i_interface],
                                                                          ln_part_tot,
                                                                          comm);
    int         **recv_stride = NULL;
    int *send_stride = PDM_array_const_int(n_gnum, 1);
    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           send_stride,
                           interface_entity2_opp_gnum[i_interface],
                           &recv_stride,
           (void ***)      &recv_interface_entity2_opp_gnum[i_interface]);
    free(send_stride);

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(recv_stride[i_part], query_itrf_n[i_interface][i_part], "recv_stride ::");

      int n_recv_tot = 0;
      for(int i = 0; i < query_itrf_n[i_interface][i_part]; ++i) {
        n_recv_tot += recv_stride[i_part][i];
      }
      n_recv_interface_entity2_opp_gnum[i_interface][i_part] = n_recv_tot;
      free(recv_stride[i_part]);
      PDM_log_trace_array_long(recv_interface_entity2_opp_gnum[i_interface][i_part], n_recv_tot, "recv_interface_entity2_opp_gnum ::");
    }
    free(recv_stride);

    PDM_block_to_part_free(btp);
  }

  for(int i_itrf = 0; i_itrf < n_interface; ++i_itrf) {
    free(interface_entity2_opp_gnum[i_itrf]);
    PDM_part_to_block_free(ptb_interface_entity2[i_itrf]);
  }
  free(interface_entity2_opp_gnum);
  free(ptb_interface_entity2);

  int          *pn_entity2_ext_opp_gnum_and_itrf = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_ext_opp_gnum_and_itrf   = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int n_recv_tot = 0;
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      n_recv_tot += n_recv_interface_entity2_opp_gnum[i_interface][i_part];
    }
    pentity2_ext_opp_gnum_and_itrf[i_part] = malloc(2 * n_recv_tot * sizeof(PDM_g_num_t));

    n_recv_tot = 0;
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      for(int i = 0; i < n_recv_interface_entity2_opp_gnum[i_interface][i_part]; ++i) {
        PDM_g_num_t gnum = PDM_ABS (recv_interface_entity2_opp_gnum[i_interface][i_part][i]);
        int         sgn  = PDM_SIGN(recv_interface_entity2_opp_gnum[i_interface][i_part][i]);
        pentity2_ext_opp_gnum_and_itrf[i_part][2*n_recv_tot  ] = gnum;
        pentity2_ext_opp_gnum_and_itrf[i_part][2*n_recv_tot+1] = -sgn * (i_interface+1); // Avec le sgn ?
        n_recv_tot++;
      }
    }

    int *order = malloc(n_recv_tot * sizeof(int));
    pn_entity2_ext_opp_gnum_and_itrf[i_part] = PDM_order_inplace_unique_and_order_long(n_recv_tot,
                                                                                   2,
                                                                                   pentity2_ext_opp_gnum_and_itrf[i_part],
                                                                                   order);
    pentity2_ext_opp_gnum_and_itrf[i_part] = realloc(pentity2_ext_opp_gnum_and_itrf[i_part], 2 * pn_entity2_ext_opp_gnum_and_itrf[i_part] * sizeof(PDM_g_num_t));


    free(order);


    PDM_log_trace_array_long(pentity2_ext_opp_gnum_and_itrf[i_part], 2 * pn_entity2_ext_opp_gnum_and_itrf[i_part], "pentity2_opp_gnum_and_itrf ::");

  }


  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      free(query_itrf_gnum   [i_interface][i_part]);
      free(recv_itrf_opp_n   [i_interface][i_part]);
      free(recv_itrf_opp_gnum[i_interface][i_part]);
      free(recv_interface_entity2_opp_gnum[i_interface][i_part]);
    }
    free(query_itrf_gnum   [i_interface]);
    free(recv_itrf_opp_n   [i_interface]);
    free(recv_itrf_opp_gnum[i_interface]);
    free(query_itrf_n      [i_interface]);
    free(recv_interface_entity2_opp_gnum  [i_interface]);
    free(n_recv_interface_entity2_opp_gnum[i_interface]);
  }
  free(recv_interface_entity2_opp_gnum  );
  free(n_recv_interface_entity2_opp_gnum);

  free(query_itrf_gnum   );
  free(query_itrf_n      );
  free(recv_itrf_opp_n   );
  free(recv_itrf_opp_gnum);




  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3,
                                                     ln_part_tot,
                                                     PDM_TRUE,
                                                     1.e-6,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_entity2, 2);


  /*
   * Prepare for unify
   */
  int         **entity2_order            = malloc( ln_part_tot * sizeof(int         *));
  PDM_g_num_t **pentity2_ln_to_gn_sorted = malloc( ln_part_tot * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    /*
     * Sort current part entity2_ln_to_gn
     */
    entity2_order           [i_part] = malloc( pn_entity2[i_part] * sizeof(int        ));
    pentity2_ln_to_gn_sorted[i_part] = malloc( pn_entity2[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2[i_part]; ++i_entity2) {
      pentity2_ln_to_gn_sorted[i_part][i_entity2] = pentity2_ln_to_gn[i_part][i_entity2];
      entity2_order           [i_part][i_entity2] = i_entity2;
    }
    PDM_sort_long(pentity2_ln_to_gn_sorted[i_part], entity2_order[i_part], pn_entity2[i_part]);
  }


  /*
   * We have two kind of extented :
   *   - From partition
   *   - From interfaces (request new gnum génération ...)
   */
  int          *pn_entity2_extented_by_interface         = malloc(ln_part_tot * sizeof(int          ));
  int          *pn_entity2_extented_by_partition         = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_extented_ln_to_gn_by_interface  = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pentity2_extented_ln_to_gn_by_partition  = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity2_extented_triplet_by_interface   = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_triplet_by_partition   = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_interface_by_interface = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_interface_by_partition = malloc(ln_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t* _pentity2_opp_gnum_and_itrf = pentity2_opp_gnum_and_itrf[i_part];

    int         *_pextract_entity1_interface = pentity1_extented_to_pentity1_interface[i_part];
    // PDM_g_num_t *_pextract_entity1_gnum      = pextract_entity1_gnum     [i_part];
    // int         *_pextract_entity1_triplet   = pextract_entity1_triplet  [i_part];

    int         *_pextract_entity1_entity2_idx     = pextract_entity1_entity2_idx    [i_part];
    PDM_g_num_t *_pextract_entity1_entity2_gnum    = pextract_entity1_entity2_gnum   [i_part];
    int         *_pextract_entity1_entity2_triplet = pextract_entity1_entity2_triplet[i_part];

    PDM_g_num_t *_pentity2_ln_to_gn_sorted   = pentity2_ln_to_gn_sorted[i_part];

    pn_entity2_extented_by_interface[i_part] = 0;
    pn_entity2_extented_by_partition[i_part] = 0;

    /* Count */
    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
      if(_pextract_entity1_interface[i_ref] != 0) {

        /* Blinder d'assert */
        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          PDM_g_num_t gnum_to_find[2] = {entity2_g_num, _pextract_entity1_interface[i_ref]};
          int pos = PDM_order_binary_search_long(gnum_to_find, _pentity2_opp_gnum_and_itrf, 2, pn_entity2_opp_gnum_and_itrf[i_part]);
          // log_trace("gnum_to_find = %i/%i --> pos = %i \n", gnum_to_find[0], gnum_to_find[1], pos);
          if(pos == -1) {
            pn_entity2_extented_by_interface[i_part]++;
          }
          // printf("Search (%i/%i) -> pos = %i \n", entity2_g_num, _pextract_entity1_interface[i_ref], pos);
        }
      } else { // Second case : Interior entity1

        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          int pos_int = PDM_binary_search_long(entity2_g_num, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);

          if(pos_int == -1) { // Not in current partition
            pn_entity2_extented_by_partition[i_part]++;
          }
        }
      }
    }

    /* Fill */
    pentity2_extented_ln_to_gn_by_interface [i_part] = malloc(2 * pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    pentity2_extented_ln_to_gn_by_partition [i_part] = malloc(    pn_entity2_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
    pentity2_extented_triplet_by_interface  [i_part] = malloc(3 * pn_entity2_extented_by_interface[i_part] * sizeof(int        ));
    pentity2_extented_triplet_by_partition  [i_part] = malloc(3 * pn_entity2_extented_by_partition[i_part] * sizeof(int        ));
    pentity2_extented_interface_by_interface[i_part] = malloc(    pn_entity2_extented_by_interface[i_part] * sizeof(int        ));
    pentity2_extented_interface_by_partition[i_part] = malloc(    pn_entity2_extented_by_partition[i_part] * sizeof(int        ));

    pn_entity2_extented_by_interface[i_part] = 0;
    pn_entity2_extented_by_partition[i_part] = 0;
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
      if(_pextract_entity1_interface[i_ref] != 0) {

        /* Blinder d'assert */
        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          PDM_g_num_t gnum_to_find[2] = {entity2_g_num, _pextract_entity1_interface[i_ref]};
          int pos = PDM_order_binary_search_long(gnum_to_find, _pentity2_opp_gnum_and_itrf, 2, pn_entity2_opp_gnum_and_itrf[i_part]);
          if(pos == -1) {
            int idx_write = pn_entity2_extented_by_interface[i_part]++;
            pentity2_extented_ln_to_gn_by_interface[i_part][2*idx_write  ] = entity2_g_num;
            pentity2_extented_ln_to_gn_by_interface[i_part][2*idx_write+1] = _pextract_entity1_interface[i_ref];

            pentity2_extented_triplet_by_interface  [i_part][3*idx_write  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
            pentity2_extented_triplet_by_interface  [i_part][3*idx_write+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
            pentity2_extented_triplet_by_interface  [i_part][3*idx_write+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
            pentity2_extented_interface_by_interface[i_part][  idx_write  ] = _pextract_entity1_interface[i_ref];
          }
        }
      } else { // Second case : Interior entity1

        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          int pos_int = PDM_binary_search_long(entity2_g_num, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);

          if(pos_int == -1) { // Not in current partition
            int idx_write = pn_entity2_extented_by_partition[i_part]++;
            pentity2_extented_ln_to_gn_by_partition [i_part][idx_write] = entity2_g_num;
            pentity2_extented_triplet_by_partition  [i_part][3*idx_write  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
            pentity2_extented_triplet_by_partition  [i_part][3*idx_write+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
            pentity2_extented_triplet_by_partition  [i_part][3*idx_write+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
            pentity2_extented_interface_by_partition[i_part][  idx_write  ] = 0; // Because interior
          }
        }
      }
    }

    PDM_log_trace_array_long(pentity2_extented_ln_to_gn_by_interface[i_part],
                             2 * pn_entity2_extented_by_interface[i_part],
                             "pentity2_extented_ln_to_gn_by_interface ::");

    PDM_gnum_set_from_parents(gen_gnum_entity2,
                              i_part,
                              pn_entity2_extented_by_interface[i_part],
                              pentity2_extented_ln_to_gn_by_interface[i_part]);

  }

  PDM_gnum_compute(gen_gnum_entity2);

  /*
   * Unify entity1_entity2
   */
  int **pextract_entity1_entity2 = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pentity2_opp_gnum_and_itrf  = pentity2_opp_gnum_and_itrf             [i_part];
    int         *_pextract_entity1_interface  = pentity1_extented_to_pentity1_interface[i_part];
    // PDM_g_num_t *_pextract_entity1_gnum      = pextract_entity1_gnum     [i_part];

    int         *_pextract_entity1_entity2_idx     = pextract_entity1_entity2_idx    [i_part];
    PDM_g_num_t *_pextract_entity1_entity2_gnum    = pextract_entity1_entity2_gnum   [i_part];
    int         *_pextract_entity1_entity2_triplet = pextract_entity1_entity2_triplet[i_part];

    int           n_entity2_opp_position     = pn_entity2_opp_gnum_and_itrf[i_part];

    PDM_g_num_t *_pentity2_ln_to_gn_sorted   = pentity2_ln_to_gn_sorted[i_part];
    int         *_entity2_order              = entity2_order           [i_part];

    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;


    pextract_entity1_entity2[i_part] = malloc( _pextract_entity1_entity2_idx[n_part1_to_part2] * sizeof(int));
    int *_pextract_entity1_entity2  = pextract_entity1_entity2 [i_part];

    PDM_g_num_t* extented_from_itrf_entity2_ln_to_gn = PDM_gnum_get(gen_gnum_entity2, i_part);

    // Erase and realloc :
    pentity2_extented_ln_to_gn_by_interface[i_part] = realloc(pentity2_extented_ln_to_gn_by_interface[i_part], pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      pentity2_extented_ln_to_gn_by_interface[i_part][i_entity2] = extented_from_itrf_entity2_ln_to_gn[i_entity2] + shift_by_domain_entity2[n_domain];
    }

    // PDM_log_trace_array_long(extented_from_itrf_entity2_ln_to_gn, pn_entity2_extented_by_interface[i_part], "extented_from_itrf_entity2_ln_to_gn ::");

    /* Sort unique the new gnum to unify */
    PDM_g_num_t *extented_from_itrf_entity2_ln_to_gn_sorted = malloc( pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      extented_from_itrf_entity2_ln_to_gn_sorted[i_entity2] = extented_from_itrf_entity2_ln_to_gn[i_entity2];
    }
    // PDM_sort_long(extented_from_itrf_entity2_ln_to_gn_sorted, extented_from_itrf_entity2_order, pn_entity2_extented_by_interface[i_part]);
    pn_entity2_extented_by_interface[i_part] = PDM_inplace_unique_long(extented_from_itrf_entity2_ln_to_gn_sorted,
                                                                       NULL,
                                                                       0,
                                                                       pn_entity2_extented_by_interface[i_part]-1);
    extented_from_itrf_entity2_ln_to_gn_sorted = realloc(extented_from_itrf_entity2_ln_to_gn_sorted, pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    int         *extented_from_itrf_entity2_order = malloc( pn_entity2_extented_by_interface[i_part] * sizeof(int        ));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      extented_from_itrf_entity2_order[i_entity2] = -1;
    }

    printf("pn_entity2_extented_by_interface[i_part] = %i \n", pn_entity2_extented_by_interface[i_part]);


    /* Sort unique the new gnum to unify */
    PDM_g_num_t *extented_from_part_entity2_ln_to_gn_sorted = malloc( pn_entity2_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
      extented_from_part_entity2_ln_to_gn_sorted[i_entity2] = pentity2_extented_ln_to_gn_by_partition[i_part][i_entity2];
      // extented_from_part_entity2_order          [i_entity2] = i_entity2;
    }
    // PDM_sort_long(extented_from_part_entity2_ln_to_gn_sorted, extented_from_part_entity2_order, pn_entity2_extented_by_partition[i_part]);

    pn_entity2_extented_by_partition[i_part] = PDM_inplace_unique_long(extented_from_part_entity2_ln_to_gn_sorted,
                                                                       NULL,
                                                                       0,
                                                                       pn_entity2_extented_by_partition[i_part]-1);
    int *extented_from_part_entity2_order = malloc( pn_entity2_extented_by_partition[i_part] * sizeof(int        ));
    extented_from_part_entity2_ln_to_gn_sorted = realloc(extented_from_part_entity2_ln_to_gn_sorted, pn_entity2_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
      extented_from_part_entity2_order[i_entity2] = -1;
    }

    printf("pn_entity2_extented_by_partition[i_part] = %i \n", pn_entity2_extented_by_partition[i_part]);


    PDM_log_trace_array_long(extented_from_part_entity2_ln_to_gn_sorted, pn_entity2_extented_by_partition[i_part], "extented_from_part_entity2_ln_to_gn_sorted ::");
    PDM_log_trace_array_long(extented_from_itrf_entity2_ln_to_gn_sorted, pn_entity2_extented_by_interface[i_part], "extented_from_itrf_entity2_ln_to_gn_sorted ::");



    /* */
    int pn_entity2_extented_by_interface2 = 0; // To read in extented_from_itrf_entity2_ln_to_gn
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {

      log_trace("_pextract_entity1_interface[i_ref] = %i \n", _pextract_entity1_interface[i_ref] );

      /*
       * First case :
       *   - entity1 is move by interface
       */
      if(_pextract_entity1_interface[i_ref] != 0) {
        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          PDM_g_num_t gnum_to_find[2] = {entity2_g_num, _pextract_entity1_interface[i_ref]};

          int pos = PDM_order_binary_search_long(gnum_to_find, _pentity2_opp_gnum_and_itrf, 2, n_entity2_opp_position);
          // log_trace("\t Search %i/%i in _pentity2_opp_gnum_and_itrf -> pos = %i \n", gnum_to_find[0], gnum_to_find[1], pos);

          if(pos == -1) {
            /*
             * Subcase :
             *   - entity2 is not in table of interface : it's a new entity2
             */
            PDM_g_num_t entity2_extented_g_num = extented_from_itrf_entity2_ln_to_gn[pn_entity2_extented_by_interface2++];
            int pos_new = PDM_binary_search_long(entity2_extented_g_num,
                                                 extented_from_itrf_entity2_ln_to_gn_sorted,
                                                 pn_entity2_extented_by_interface[i_part]);

            assert(pos_new != -1);
            int shift = pn_entity2[i_part] + pn_entity2_extented_by_partition[i_part];

            if(extented_from_itrf_entity2_order[pos_new] == -1) {
              extented_from_itrf_entity2_order[pos_new] = 1;

              pentity2_extented_ln_to_gn_by_interface [i_part][  pos_new  ] = entity2_extented_g_num + shift_by_domain_entity2[n_domain];

              pentity2_extented_triplet_by_interface  [i_part][3*pos_new  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
              pentity2_extented_triplet_by_interface  [i_part][3*pos_new+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
              pentity2_extented_triplet_by_interface  [i_part][3*pos_new+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
              pentity2_extented_interface_by_interface[i_part][  pos_new  ] = _pextract_entity1_interface[i_ref];

            }

            // _pextract_entity1_entity2[idx_entity1] = ( shift + extented_from_itrf_entity2_order[pos_new] + 1); // ATTENTION SIGN
            _pextract_entity1_entity2[idx_entity1] = ( shift + pos_new + 1); // ATTENTION SIGN

            log_trace("\t Translate cas 1 : %i  ---> (idx_entity1 = %i) --> %i \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1]);

          } else {
            /*
             * Subcase :
             *   - entity2 is in table of interface
             */
            int pos2      = entity2_opp_position[i_part][pos ];
            _pextract_entity1_entity2[idx_entity1] = ( pos2 + 1);

            // log_trace("\t Translate cas 2 : %i  ---> (idx_entity1 = %i) --> %i (pos=%i/pos2=%i) \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1], pos, pos2);


          }
        }
      } else { // Second case : Interior entity1

        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          int pos_int = PDM_binary_search_long(entity2_g_num, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);

          if(pos_int != -1) {                   // entity2 is already in current partition
            assert(pos_int != -1);
            int orig_pos = _entity2_order[pos_int];
            _pextract_entity1_entity2[idx_entity1] = ( orig_pos + 1);

            log_trace("\t Translate cas 3 : %i  ---> (idx_entity1 = %i) --> %i \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1]);

          } else {                              // Not in current partition
            int pos_ext = PDM_binary_search_long(entity2_g_num,
                                                 extented_from_part_entity2_ln_to_gn_sorted,
                                                 pn_entity2_extented_by_partition[i_part]);
            assert(pos_ext != -1);

            if(extented_from_part_entity2_order[pos_ext] == -1) {
              extented_from_part_entity2_order[pos_ext] = 1;
              pentity2_extented_ln_to_gn_by_partition [i_part][  pos_ext  ] = entity2_g_num;

              pentity2_extented_triplet_by_partition  [i_part][3*pos_ext  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
              pentity2_extented_triplet_by_partition  [i_part][3*pos_ext+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
              pentity2_extented_triplet_by_partition  [i_part][3*pos_ext+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
              pentity2_extented_interface_by_partition[i_part][  pos_ext  ] = 0;

            }

            // _pextract_entity1_entity2[idx_entity1] = ( pn_entity2[i_part] + extented_from_part_entity2_order[pos_ext] + 1);
            _pextract_entity1_entity2[idx_entity1] = ( pn_entity2[i_part] + pos_ext + 1);

            log_trace("\t Translate cas 4 : %i  ---> (idx_entity1 = %i) --> %i \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1]);

          }
        }
      }
    }

    // assert(pn_entity2_extented_by_interface2 == pn_entity2_extented_by_interface[i_part]);

    if(1 == 1) {
      PDM_log_trace_connectivity_int(_pextract_entity1_entity2_idx,
                                     _pextract_entity1_entity2,
                                     n_part1_to_part2,
                                     "pextract_entity1_entity2 ::");
    }

    free(extented_from_itrf_entity2_order          );
    free(extented_from_itrf_entity2_ln_to_gn_sorted);
    free(extented_from_part_entity2_order          );
    free(extented_from_part_entity2_ln_to_gn_sorted);

  }


  PDM_gnum_free(gen_gnum_entity2);

  /*
   * Reconstruction
   */
  int          *pn_entity2_extented              = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_extented_ln_to_gn       = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity2_extented_triplet        = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_interface      = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_to_entity2_idx = malloc(ln_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    pn_entity2_extented[i_part] = pn_entity2_extented_by_partition[i_part] + pn_entity2_extented_by_interface[i_part];
    pentity2_extented_ln_to_gn      [i_part] = malloc(    pn_entity2_extented[i_part]      * sizeof(PDM_g_num_t));
    pentity2_extented_triplet       [i_part] = malloc(3 * pn_entity2_extented[i_part]      * sizeof(int        ));
    pentity2_extented_interface     [i_part] = malloc(    pn_entity2_extented[i_part]      * sizeof(int        ));
    pentity2_extented_to_entity2_idx[i_part] = malloc( (  pn_entity2_extented[i_part] + 1) * sizeof(int        ));

    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
      pentity2_extented_ln_to_gn [i_part][  i_entity2  ] = pentity2_extented_ln_to_gn_by_partition [i_part][  i_entity2  ];
      pentity2_extented_triplet  [i_part][3*i_entity2  ] = pentity2_extented_triplet_by_partition  [i_part][3*i_entity2  ];
      pentity2_extented_triplet  [i_part][3*i_entity2+1] = pentity2_extented_triplet_by_partition  [i_part][3*i_entity2+1];
      pentity2_extented_triplet  [i_part][3*i_entity2+2] = pentity2_extented_triplet_by_partition  [i_part][3*i_entity2+2];
      pentity2_extented_interface[i_part][  i_entity2  ] = pentity2_extented_interface_by_partition[i_part][  i_entity2  ];
    }

    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      int idx_write = pn_entity2_extented_by_partition[i_part]+i_entity2;
      pentity2_extented_ln_to_gn [i_part][  idx_write  ] = pentity2_extented_ln_to_gn_by_interface [i_part][  i_entity2  ];
      pentity2_extented_triplet  [i_part][3*idx_write  ] = pentity2_extented_triplet_by_interface  [i_part][3*i_entity2  ];
      pentity2_extented_triplet  [i_part][3*idx_write+1] = pentity2_extented_triplet_by_interface  [i_part][3*i_entity2+1];
      pentity2_extented_triplet  [i_part][3*idx_write+2] = pentity2_extented_triplet_by_interface  [i_part][3*i_entity2+2];
      pentity2_extented_interface[i_part][  idx_write  ] = pentity2_extented_interface_by_interface[i_part][  i_entity2  ];
    }

    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented[i_part]+1; ++i_entity2) {
      pentity2_extented_to_entity2_idx[i_part][i_entity2] = 3 * i_entity2;
    }

    if(1 == 1) {
      PDM_log_trace_array_int (pentity2_extented_to_entity2_idx[i_part],     pn_entity2_extented[i_part]+1, "pentity2_extented_to_entity2_idx ::" );
      PDM_log_trace_array_int (pentity2_extented_triplet       [i_part], 3 * pn_entity2_extented[i_part]  , "pentity2_extented_triplet    ::" );
      PDM_log_trace_array_int (pentity2_extented_interface     [i_part],     pn_entity2_extented[i_part]  , "pentity2_extented_interface  ::" );
      PDM_log_trace_array_long(pentity2_extented_ln_to_gn      [i_part],     pn_entity2_extented[i_part]  , "pentity2_extented_ln_to_gn   ::" );
    }

  }

  *pn_entity2_extented_out                     = pn_entity2_extented;
  *pentity2_extented_ln_to_gn_out              = pentity2_extented_ln_to_gn;
  *pextented_entity1_entity2_idx_out           = pextract_entity1_entity2_idx;
  *pextented_entity1_entity2_out               = pextract_entity1_entity2;
  *pentity2_extented_to_pentity2_idx_out       = pentity2_extented_to_entity2_idx;
  *pentity2_extented_to_pentity2_triplet_out   = pentity2_extented_triplet;
  *pentity2_extented_to_pentity2_interface_out = pentity2_extented_interface;



  // for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
  //   free(pextract_entity1_entity2[i_part]);
  // }
  // free(pextract_entity1_entity2);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_entity1_entity2_n      [i_part]);
    free(pextract_entity1_entity2_gnum   [i_part]);
    free(pextract_entity1_entity2_triplet[i_part]);
    // free(pextract_entity1_entity2_idx [i_part]);
  }
  free(pextract_entity1_entity2_n      );
  free(pextract_entity1_entity2_gnum   );
  free(pextract_entity1_entity2_triplet);
  // free(pextract_entity1_entity2_idx );


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity2_extented_ln_to_gn_by_interface [i_part]);
    free(pentity2_extented_ln_to_gn_by_partition [i_part]);
    free(pentity2_extented_triplet_by_interface  [i_part]);
    free(pentity2_extented_triplet_by_partition  [i_part]);
    free(pentity2_extented_interface_by_interface[i_part]);
    free(pentity2_extented_interface_by_partition[i_part]);
  }
  free(pentity2_extented_ln_to_gn_by_interface );
  free(pentity2_extented_ln_to_gn_by_partition );
  free(pentity2_extented_triplet_by_interface  );
  free(pentity2_extented_triplet_by_partition  );
  free(pentity2_extented_interface_by_interface);
  free(pentity2_extented_interface_by_partition);

  free(pn_entity2_extented_by_interface );
  free(pn_entity2_extented_by_partition );

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(entity2_order           [i_part]);
    free(pentity2_ln_to_gn_sorted[i_part]);
  }
  free(entity2_order           );
  free(pentity2_ln_to_gn_sorted);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity2_opp_gnum_and_itrf[i_part]);
    free(pentity2_ext_opp_gnum_and_itrf[i_part]);
    free(entity2_opp_position      [i_part]);
  }
  free(entity2_opp_position        );
  free(pentity2_opp_gnum_and_itrf  );
  free(pentity2_ext_opp_gnum_and_itrf  );
  free(pn_entity2_opp_gnum_and_itrf);
  free(pn_entity2_ext_opp_gnum_and_itrf);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity2_num             [i_part]);
    free(pentity2_opp_location_idx[i_part]);
    free(pentity2_opp_location    [i_part]);
    free(pentity2_opp_interface   [i_part]);
    free(pentity2_opp_sens        [i_part]);
    free(pentity2_opp_gnum        [i_part]);
  }
  free(pn_entity2_num            );
  free(pentity2_num              );
  free(pentity2_opp_location_idx );
  free(pentity2_opp_location     );
  free(pentity2_opp_interface_idx);
  free(pentity2_opp_interface    );
  free(pentity2_opp_sens         );
  free(pentity2_opp_gnum         );

  free(pn_entity1          );
  free(pentity1_ln_to_gn   );
  free(pentity1_entity2_idx);
  free(pentity1_entity2    );
  free(pn_entity2          );
  free(pentity2_ln_to_gn   );

}

static void
PDM_part_extension_compute2
(
        PDM_part_extension_t *part_ext,
  const int                   dim
)
{
  // TODO : mv dim in create but break API

  /* Manage shift */
  // _offset_parts_by_domain(part_ext, 1);

  /* Manage dim */
  if(dim == 3) {
    _part_extension_3d(part_ext);
  } else if(dim == 2) {
    _part_extension_2d(part_ext);
  } else if(dim == 1) {
    _part_extension_1d(part_ext);
  } else  {
    PDM_error(__FILE__, __LINE__, 0, "Wrong dim size in PDM_part_extension_compute2 : %d ( Should be >=1 )\n", (int) dim);
  }


  /* Manage loop for depth AND multiple transformation */



  /* Manage unshift */
  // _offset_parts_by_domain(part_ext, -1);
  // _offset_results_by_domain(part_ext);

}



#ifdef __cplusplus
}
#endif /* __cplusplus */

