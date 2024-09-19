#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_dcube_gen.h"
#include "pdm_distrib.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

#include "pdm_domain_interface.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_exchange_point_list
(
 int            n_group_join,
 int           *group_join_to_join_opp,
 int           *dface_join_idx,
 PDM_g_num_t   *dface_join,
 PDM_g_num_t  **dface_join_opp,
 PDM_MPI_Comm   comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * We have for all extraction domain, we need to exchange id
   */
  PDM_g_num_t **distrib_join   = malloc(n_group_join * sizeof(PDM_g_num_t *));
  // PDM_g_num_t **dface_join_opp = malloc(n_group_join * sizeof(PDM_g_num_t *));
  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    distrib_join[i_group_join] = PDM_compute_entity_distribution(comm, dn_face_join);
  }

  *dface_join_opp = malloc(dface_join_idx[n_group_join] * sizeof(PDM_g_num_t));
  PDM_g_num_t* _dface_join_opp = *dface_join_opp;

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {

    int i_group_join_opp = group_join_to_join_opp[i_group_join];
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    PDM_g_num_t* distrib_join_cur = distrib_join[i_group_join    ];
    PDM_g_num_t* distrib_join_opp = distrib_join[i_group_join_opp];

    PDM_g_num_t *join_ln_to_gn = malloc(dn_face_join * sizeof(PDM_g_num_t));
    for(int i = 0; i < dn_face_join; ++i) {
      join_ln_to_gn[i] = distrib_join_cur[i_rank] + i + 1;
    }

    /*
     * Exchange
     */
    PDM_g_num_t *blk_dface_join_cur = &dface_join[dface_join_idx[i_group_join    ]];
    PDM_g_num_t *blk_dface_join_opp = &dface_join[dface_join_idx[i_group_join_opp]];

    PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_join_opp,
                                 (const PDM_g_num_t **) &join_ln_to_gn,
                                                        &dn_face_join,
                                                        1,
                                                        comm);

    int cst_stride = 1;
    PDM_g_num_t* sub_dface_join_opp = &_dface_join_opp[dface_join_idx[i_group_join]];
    PDM_block_to_part_exch_in_place(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &cst_stride,
                  (void *) blk_dface_join_opp,
                           NULL,
               (void ** ) &sub_dface_join_opp);

    if(0 == 1) {
      PDM_log_trace_array_long(blk_dface_join_cur, dn_face_join, "dface_join_cur :: ");
      PDM_log_trace_array_long(sub_dface_join_opp, dn_face_join, "dface_join_opp :: ");
    }

    PDM_block_to_part_free(btp);
    free(join_ln_to_gn);
  }

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    free(distrib_join[i_group_join]);
  }
  free(distrib_join);
}



/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_vtx_seg,
           double        *length,
           int           *verbose)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_domain    = 2;
  int                verbose   = 0;


  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &verbose);

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Alloc for distributed mesh */
  int *dn_cell       = (int *) malloc(n_domain * sizeof(int));
  int *dn_face       = (int *) malloc(n_domain * sizeof(int));
  int *dn_vtx        = (int *) malloc(n_domain * sizeof(int));
  int *n_face_group  = (int *) malloc(n_domain * sizeof(int));
  int *dface_vtx_s   = (int *) malloc(n_domain * sizeof(int));
  int *dface_group_s = (int *) malloc(n_domain * sizeof(int));

  PDM_g_num_t  **dface_cell      = (PDM_g_num_t **) malloc(n_domain * sizeof(PDM_g_num_t *));
  int          **dface_vtx_idx   = (int         **) malloc(n_domain * sizeof(int         *));
  PDM_g_num_t  **dface_vtx       = (PDM_g_num_t **) malloc(n_domain * sizeof(PDM_g_num_t *));
  double       **dvtx_coord      = (double      **) malloc(n_domain * sizeof(double      *));
  int          **dface_group_idx = (int         **) malloc(n_domain * sizeof(int         *));
  PDM_g_num_t  **dface_group     = (PDM_g_num_t **) malloc(n_domain * sizeof(PDM_g_num_t *));


  int          **dface_bnd_idx   = (int         **) malloc(n_domain * sizeof(int         *));
  PDM_g_num_t  **dface_bnd       = (PDM_g_num_t **) malloc(n_domain * sizeof(PDM_g_num_t *));

  int n_group_join = 2*(n_domain-1);
  int          *dface_join_idx  = (int         *) malloc((n_group_join + 1) * sizeof(int        ));
  PDM_g_num_t  *dface_join      = NULL; // A allouer propremet

  int *group_join_to_domain_cur = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_domain_opp = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_join_opp = (int *) malloc( n_group_join * sizeof(int));


  int tmp_i_domain = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    if (i_join % 2 == 0) {
      group_join_to_join_opp[i_join] = i_join + 1;
      group_join_to_domain_opp[i_join] = tmp_i_domain + 1;
    } else {
      group_join_to_join_opp[i_join] = i_join - 1;
      group_join_to_domain_opp[i_join] = tmp_i_domain++;
    }
  }


  for(int i_group_join = 0; i_group_join < n_group_join+1; ++i_group_join) {
    dface_join_idx[i_group_join] = 0;
  }

  PDM_dcube_t **dcube = (PDM_dcube_t **) malloc(n_domain * sizeof(PDM_dcube_t *));
  int tmp_i_group_join = 0;
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {

    dcube[i_domain] = PDM_dcube_gen_init(comm, n_vtx_seg, length, i_domain, 0., 0., PDM_OWNERSHIP_KEEP);
    PDM_dcube_gen_dim_get(dcube         [i_domain],
                          &n_face_group [i_domain],
                          &dn_cell      [i_domain],
                          &dn_face      [i_domain],
                          &dn_vtx       [i_domain],
                          &dface_vtx_s  [i_domain],
                          &dface_group_s[i_domain]);

    PDM_dcube_gen_data_get(dcube          [i_domain],
                          &dface_cell     [i_domain],
                          &dface_vtx_idx  [i_domain],
                          &dface_vtx      [i_domain],
                          &dvtx_coord     [i_domain],
                          &dface_group_idx[i_domain],
                          &dface_group    [i_domain]);

    /*
     * Les faces groups du dcube sont : zmin, zmax, xmin, xmax, ymin, ymax
     * Il faut les séparer en faces de bords et faces raccord, sachant que
     * les domains sont alignées selon X
     */
    int n_bnd = 4;
    int n_jn  = 2;
    if (i_domain == 0){
      n_bnd++;
      n_jn-- ;
    }
    if (i_domain == n_domain-1){
      n_bnd++;
      n_jn-- ;
    }
    // Join numbering (left to right, increasing i_domain)

    dface_bnd_idx [i_domain] = (int *) malloc((n_bnd        + 1) * sizeof(int));

    // First pass to count and allocate
    int i_bnd = 1;
    int i_jn  = tmp_i_group_join+1;
    dface_bnd_idx[i_domain][0]  = 0;
    // dface_join_idx[0] = 0;
    for (int igroup = 0; igroup < n_face_group[i_domain]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_domain == 0)) && (igroup != 3 || (igroup == 3 && i_domain == n_domain-1));
      int group_size = dface_group_idx[i_domain][igroup+1] - dface_group_idx[i_domain][igroup];
      if (copy_to_bnd) { //Its a boundary
        dface_bnd_idx[i_domain][i_bnd++] = group_size;
      } else { //Its a join
        group_join_to_domain_cur[i_jn-1] = i_domain;
        dface_join_idx[i_jn++] = group_size;
      }
    }
    for (int i = 0; i < n_bnd; i++) {
      dface_bnd_idx[i_domain][i+1] = dface_bnd_idx[i_domain][i+1] + dface_bnd_idx[i_domain][i];
    }


    for (int i = tmp_i_group_join; i < i_jn-1; i++) {
      dface_join_idx[i+1] = dface_join_idx[i+1] + dface_join_idx[i];
    }

    /* A bit stupid but complicated to made it in ohter way for a clear test */
    dface_join = realloc(dface_join, dface_join_idx[i_jn-1] * sizeof(PDM_g_num_t));

    // Second pass to copy
    dface_bnd [i_domain] = (PDM_g_num_t *) malloc(dface_bnd_idx [i_domain][n_bnd        ] * sizeof(PDM_g_num_t));
    i_bnd = 0;
    i_jn  = tmp_i_group_join;
    for (int igroup = 0; igroup < n_face_group[i_domain]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_domain == 0)) && (igroup != 3 || (igroup == 3 && i_domain == n_domain-1));
      if (copy_to_bnd){ //Its a boundary
        for (int i = dface_group_idx[i_domain][igroup]; i < dface_group_idx[i_domain][igroup+1]; i++) {
          dface_bnd[i_domain][i_bnd++] = dface_group[i_domain][i];
        }
      } else { //Its a join
        int k = 0;
        for (int i = dface_group_idx[i_domain][igroup]; i < dface_group_idx[i_domain][igroup+1]; i++) {
          dface_join[dface_join_idx[i_jn]+k++] = dface_group[i_domain][i];
        }
        i_jn++;
      }
    }

    /*
     *  Go to nexts join
     */
    tmp_i_group_join += n_jn;
  }

  /*
   * Setup dface_join_opp + group_id in current layout
   */
  PDM_g_num_t *dface_join_opp = NULL;
  _exchange_point_list(n_group_join,
                       group_join_to_join_opp,
                       dface_join_idx,
                       dface_join,
                       &dface_join_opp,
                       comm);

  /*log_trace("Global join data (%d)\n", n_group_join);*/
  /*PDM_log_trace_array_int(group_join_to_domain_cur, n_group_join, "group_join_to_domain_cur :: ");*/
  /*PDM_log_trace_array_int(group_join_to_join_opp, n_group_join, "group_join_to_join_opp :: ");*/
  /*PDM_log_trace_array_int(group_join_to_domain_opp, n_group_join, "group_join_to_domain_opp :: ");*/

  /*PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "dface_join_idx :: ");*/
  /*PDM_log_trace_array_long(dface_join    , dface_join_idx[n_group_join], "dface_join     :: ");*/
  /*PDM_log_trace_array_long(dface_join_opp, dface_join_idx[n_group_join], "dface_join_opp :: ");*/

  // Convert for new version
  int n_interface = n_group_join / 2;
  int          *interface_dn_f  = malloc(n_interface * sizeof(int));
  PDM_g_num_t **interface_ids_f = malloc(n_interface * sizeof(PDM_g_num_t*));
  int         **interface_dom_f = malloc(n_interface * sizeof(int*));

  int i_interface = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    if (i_join <= i_join_opp) {
      int n_face = dface_join_idx[i_join+1] - dface_join_idx[i_join];
      interface_dn_f [i_interface] = n_face;
      interface_ids_f[i_interface] = malloc(2*n_face*sizeof(PDM_g_num_t));
      interface_dom_f[i_interface] = malloc(2*n_face*sizeof(int));
      int idx = 0;
      for (int i_face = dface_join_idx[i_join]; i_face < dface_join_idx[i_join+1]; i_face++) {
        interface_ids_f[i_interface][idx] = dface_join[i_face];
        interface_dom_f[i_interface][idx++] = group_join_to_domain_cur[i_join];

        interface_ids_f[i_interface][idx] = dface_join_opp[i_face];
        interface_dom_f[i_interface][idx++] = group_join_to_domain_opp[i_join];
      }
      assert (idx == 2*n_face);
      i_interface++;
    }
  }
  assert (i_interface == n_interface);

  if (verbose) {
    log_trace("Number of interfaces : %d\n", n_interface);
    PDM_log_trace_array_int(interface_dn_f, n_interface, "interface_dn_f ::");
    for (i_interface = 0; i_interface < n_interface; i_interface ++) {
      log_trace("Interface %d\n", i_interface);
      PDM_log_trace_array_long(interface_ids_f[i_interface], 2*interface_dn_f[i_interface], "  face ids ::");
      PDM_log_trace_array_int (interface_dom_f[i_interface], 2*interface_dn_f[i_interface], "  face dom ::");
    }
  }
  
  // New version begins
  PDM_domain_interface_t *dom_intrf = PDM_domain_interface_create(n_interface, n_domain,
      PDM_DOMAIN_INTERFACE_MULT_YES, PDM_OWNERSHIP_KEEP, comm);
  PDM_domain_interface_set(dom_intrf, PDM_BOUND_TYPE_FACE, interface_dn_f, interface_ids_f, interface_dom_f);

  int          *interface_dn_v  = NULL;
  PDM_g_num_t **interface_ids_v = NULL;
  int         **interface_dom_v = NULL;
  PDM_domain_interface_translate_face2vtx(dom_intrf, dn_vtx, dn_face, dface_vtx_idx, dface_vtx);
  PDM_domain_interface_get(dom_intrf, PDM_BOUND_TYPE_VTX, &interface_dn_v, &interface_ids_v, &interface_dom_v);

  if (verbose) {
    PDM_log_trace_array_int(interface_dn_v, n_interface, "interface_dn_v ::");
    for (i_interface = 0; i_interface < n_interface; i_interface ++) {
      log_trace("Interface %d\n", i_interface);
      PDM_log_trace_array_long(interface_ids_v[i_interface], 2*interface_dn_v[i_interface], "  vtx ids ::");
      PDM_log_trace_array_int (interface_dom_v[i_interface], 2*interface_dn_v[i_interface], "  vtx dom ::");
    }
  }


  PDM_domain_interface_free(dom_intrf);

  /* Free memory */
  for (int i_domain = 0; i_domain < n_domain; i_domain++) {
    free(dface_bnd_idx [i_domain]);
    free(dface_bnd     [i_domain]);
    PDM_dcube_gen_free(dcube[i_domain]);
  }
  free(dcube);
  free(dn_cell);
  free(dn_face);
  free(dn_vtx);
  free(n_face_group);
  free(dface_group_s);
  free(dface_vtx_s);
  free(dface_cell);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(dvtx_coord);
  free(dface_group_idx);
  free(dface_group);
  free(dface_bnd_idx);
  free(dface_bnd);
  free(dface_join_idx);
  free(dface_join);
  free(dface_join_opp);

  free(group_join_to_domain_cur);
  free(group_join_to_domain_opp);
  free(group_join_to_join_opp);

  for (i_interface = 0; i_interface < n_interface; i_interface++)
  {
    free(interface_ids_f[i_interface]);
    free(interface_dom_f[i_interface]);
  }
  free(interface_dn_f );
  free(interface_ids_f);
  free(interface_dom_f);

  PDM_MPI_Finalize();

  return 0;
}
