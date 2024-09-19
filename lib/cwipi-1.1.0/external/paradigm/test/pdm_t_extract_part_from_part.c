#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dcube_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_gnum.h"
#include "pdm_part_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_dmesh.h"
#include "pdm_unique.h"
#include "pdm_part_geom.h"
#include "pdm_logging.h"
#include "pdm_priv.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

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
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
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
 * \param [inout]   part_method Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *part_method)
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
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
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

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method);

  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell = NULL;
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t  *dface_vtx = NULL;
  double       *dvtx_coord = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group = NULL;
  int           dface_vtxL;
  int           dFaceGroupL;

  PDM_dcube_t* dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          0.,
                                          0.,
                                          0.,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                        &n_face_group,
                        &dn_cell,
                        &dn_face,
                        &dn_vtx,
                        &dface_vtxL,
                        &dFaceGroupL);

  PDM_dcube_gen_data_get(dcube,
                         &dface_cell,
                         &dface_vtx_idx,
                         &dface_vtx,
                         &dvtx_coord,
                         &dface_group_idx,
                         &dface_group);

  /*
   * Create dmesh
   */
  PDM_dmesh_t* dm = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     dn_cell,
                                     dn_face,
                                     0, // dn_edge
                                     dn_vtx,
                                     comm);

  PDM_dmesh_vtx_coord_set(dm,
                          dvtx_coord,
                          PDM_OWNERSHIP_USER);


  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             dface_vtx,
                             dface_vtx_idx,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_FACE_CELL,
                             dface_cell,
                             NULL,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_bound_set(dm,
                      PDM_BOUND_TYPE_FACE,
                      n_face_group,
                      dface_group,
                      dface_group_idx,
                      PDM_OWNERSHIP_USER);

  /*
   * Partitionnement
   */
  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart, -1, "PDM_PART_RENUM_CELL_NONE", NULL, "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_dmesh_set(mpart, 0, dm);
  PDM_multipart_compute(mpart);

  /*
   * Get the partition domain
   */
  int i_domain = 0;

  double      **cell_center             = (double      **) malloc( n_part_domains * sizeof(double      *));
  PDM_g_num_t **selected_g_num          = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pcell_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pface_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pvtx_ln_to_gn           = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  int         **selected_g_num_idx      = (int         **) malloc( n_part_domains * sizeof(int         *));
  int          *pn_cell                 = (int          *) malloc( n_part_domains * sizeof(int          ));
  int          *pn_face                 = (int          *) malloc( n_part_domains * sizeof(int          ));
  int          *pn_vtx                  = (int          *) malloc( n_part_domains * sizeof(int          ));
  int          *pn_select_cell          = (int          *) malloc( n_part_domains * sizeof(int          ));
  double      **weight                  = (double      **) malloc( n_part_domains * sizeof(double      *));
  int         **pcell_face              = (int         **) malloc( n_part_domains * sizeof(int         *));
  int         **pcell_face_idx          = (int         **) malloc( n_part_domains * sizeof(int         *));
  int         **pface_vtx               = (int         **) malloc( n_part_domains * sizeof(int         *));
  int         **pface_vtx_idx           = (int         **) malloc( n_part_domains * sizeof(int         *));
  double      **pvtx_coord              = (double      **) malloc( n_part_domains * sizeof(double      *));
  double      **tmp_extract_cell_center = (double      **) malloc( n_part_domains * sizeof(double      *));

  PDM_gen_gnum_t* gnum_extract = PDM_gnum_create(3,
                                                 n_part_domains,
                                                 PDM_FALSE,
                                                 1.e-6,
                                                 comm,
                                                 PDM_OWNERSHIP_KEEP);

  for (int i_part = 0; i_part < n_part_domains; i_part++){

    PDM_g_num_t* cell_ln_to_gn = NULL;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    i_domain,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    &cell_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    int *cell_face     = NULL;
    int *cell_face_idx = NULL;
    int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                     i_domain,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                     &cell_face_idx,
                                                     &cell_face,
                                                     PDM_OWNERSHIP_KEEP);

    int *face_vtx     = NULL;
    int *face_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        i_domain,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &face_vtx_idx,
                                        &face_vtx,
                                        PDM_OWNERSHIP_KEEP);

    PDM_g_num_t* face_ln_to_gn = NULL;
    int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                 i_domain,
                                                 i_part,
                                                 PDM_MESH_ENTITY_FACE,
                                                 &face_ln_to_gn,
                                                 PDM_OWNERSHIP_KEEP);

    PDM_g_num_t* vtx_ln_to_gn = NULL;
    int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_domain,
                                                i_part,
                                                PDM_MESH_ENTITY_VTX,
                                                &vtx_ln_to_gn,
                                                PDM_OWNERSHIP_KEEP);

    double *vtx = NULL;
    PDM_multipart_part_vtx_coord_get(mpart,
                                     i_domain,
                                     i_part,
                                     &vtx,
                                     PDM_OWNERSHIP_KEEP);

    pn_cell       [i_part] = n_cell;
    pcell_ln_to_gn[i_part] = cell_ln_to_gn;
    pface_ln_to_gn[i_part] = face_ln_to_gn;
    pvtx_ln_to_gn [i_part] = vtx_ln_to_gn;
    pcell_face    [i_part] = cell_face;
    pcell_face_idx[i_part] = cell_face_idx;
    pn_face       [i_part] = n_face;
    pn_vtx        [i_part] = n_vtx;

    pface_vtx    [i_part] = face_vtx;
    pface_vtx_idx[i_part] = face_vtx_idx;
    pvtx_coord   [i_part] = vtx;


    /*
     * Compute center-cell and extract cells corresponding to criteria
     */
    double *face_center         = (double *) malloc( 3 * n_face * sizeof(double));

    for(int i_face = 0; i_face < n_face; ++i_face) {
      face_center[3*i_face  ] = 0.;
      face_center[3*i_face+1] = 0.;
      face_center[3*i_face+2] = 0.;
      int n_vtx_on_face = face_vtx_idx[i_face+1] - face_vtx_idx[i_face];
      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = PDM_ABS(face_vtx[idx_vtx])-1;
        face_center[3*i_face  ] += vtx[3*i_vtx  ];
        face_center[3*i_face+1] += vtx[3*i_vtx+1];
        face_center[3*i_face+2] += vtx[3*i_vtx+2];
      }
      face_center[3*i_face  ] = face_center[3*i_face  ] / n_vtx_on_face;
      face_center[3*i_face+1] = face_center[3*i_face+1] / n_vtx_on_face;
      face_center[3*i_face+2] = face_center[3*i_face+2] / n_vtx_on_face;
    }

    cell_center[i_part] = (double *) malloc( 3 * n_cell * sizeof(double));

    for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

      cell_center[i_part][3*i_cell  ] = 0.;
      cell_center[i_part][3*i_cell+1] = 0.;
      cell_center[i_part][3*i_cell+2] = 0.;

      int n_face_on_cell = cell_face_idx[i_cell+1] - cell_face_idx[i_cell];

      for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {

        int i_face = PDM_ABS(cell_face[idx_face])-1;
        cell_center[i_part][3*i_cell  ] += face_center[3*i_face  ];
        cell_center[i_part][3*i_cell+1] += face_center[3*i_face+1];
        cell_center[i_part][3*i_cell+2] += face_center[3*i_face+2];
      }
      cell_center[i_part][3*i_cell  ] = cell_center[i_part][3*i_cell  ] / n_face_on_cell;
      cell_center[i_part][3*i_cell+1] = cell_center[i_part][3*i_cell+1] / n_face_on_cell;
      cell_center[i_part][3*i_cell+2] = cell_center[i_part][3*i_cell+2] / n_face_on_cell;
    }

    free(face_center);

    selected_g_num         [i_part] = (PDM_g_num_t *) malloc(  n_cell          * sizeof(PDM_g_num_t));
    selected_g_num_idx     [i_part] = (int         *) malloc( (n_cell + 1)     * sizeof(int        ));
    tmp_extract_cell_center[i_part] = (double      *) malloc(  3 * n_cell      * sizeof(double     ));

    /*
     * Sub-part
     */

    double bbox[6];
    bbox[0] = 0.3;
    bbox[1] = 0.3;
    bbox[2] = 0.35;
    bbox[3] = 0.7;
    bbox[4] = 0.7;
    bbox[5] = 0.65;

    int n_select_cell = 0;
    selected_g_num_idx[i_part][0] = 0;
    for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

      int inside = 1;
      for(int i = 0; i < 3; ++i) {
        if (cell_center[i_part][3*i_cell+i] > bbox[i+3] || cell_center[i_part][3*i_cell+i] < bbox[i]) {
          inside = 0;
        }
      }
      if(inside == 1) {
        selected_g_num         [i_part][n_select_cell]     = cell_ln_to_gn[i_cell];
        tmp_extract_cell_center[i_part][3*n_select_cell  ] = cell_center[i_part][3*i_cell  ];
        tmp_extract_cell_center[i_part][3*n_select_cell+1] = cell_center[i_part][3*i_cell+1];
        tmp_extract_cell_center[i_part][3*n_select_cell+2] = cell_center[i_part][3*i_cell+2];
        n_select_cell++;

      }
      selected_g_num_idx[i_part][i_cell+1] = selected_g_num_idx[i_part][i_cell] + inside;
    }

    selected_g_num         [i_part] = realloc(selected_g_num[i_part], n_select_cell * sizeof(PDM_g_num_t));
    pn_select_cell         [i_part] = n_select_cell;
    tmp_extract_cell_center[i_part] = realloc(tmp_extract_cell_center[i_part], 3 * n_select_cell * sizeof(double));

    weight        [i_part] = malloc(n_select_cell * sizeof(double));
    for(int i = 0; i < n_select_cell; ++i) {
      weight        [i_part][i] = 1.;
    }

    // PDM_log_trace_array_long(selected_g_num    [i_part], n_select_cell, "selected_g_num     : ");
    // PDM_log_trace_array_int (selected_g_num_idx[i_part], n_cell+1     , "selected_g_num_idx : ");
    // PDM_log_trace_array_long(cell_ln_to_gn             , n_cell       , "cell_ln_to_gn      : ");


  }

  /*
   *  We know the extraction of cells required : now we use part_to_part to resetup a coherent part :
   *           - Redistribute over all process ( And order with scoth or something else)
   *           - Keep parent_cell_g_num
   *           - Rebuild by descending connectivity a coherent partition (cell -> face -> vtx + face_group)
   */

  PDM_g_num_t **child_selected_g_num = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  for (int i_part = 0; i_part < n_part_domains; i_part++){
    PDM_gnum_set_from_coords(gnum_extract, i_part, pn_select_cell[i_part], tmp_extract_cell_center[i_part], NULL);
  }

  PDM_gnum_compute(gnum_extract);

  for (int i_part = 0; i_part < n_part_domains; i_part++){
    child_selected_g_num[i_part] = PDM_gnum_get(gnum_extract, i_part);
    // PDM_log_trace_array_long(child_selected_g_num[i_part], pn_select_cell[i_part], "child_selected_g_num : ");
  }

  // Ameliorer par himbert ordering
  // Sinon on doit calculer un graphe distribué
  // Donc cell_cell + part_to_part + on enleve les cellules hors de la numérotation ?

  /*
   *  Remake equilibrate block -> Block is not partial
   */
  PDM_part_to_block_t *ptb_equi = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           child_selected_g_num,
                                                           weight,
                                                           pn_select_cell,
                                                           n_part_domains,
                                                           comm);

  int dn_cell_equi = PDM_part_to_block_n_elt_block_get (ptb_equi);
  PDM_g_num_t *dextract_gnum = PDM_part_to_block_block_gnum_get(ptb_equi);

  // /*
  //  * Exchange cell center to reorder by hilbert
  //  */
  // int    *blk_center_n    = NULL;
  // double *blk_cell_center = NULL;

  // int** cell_center_n = malloc(n_part_domains * sizeof(int * ));
  // for(int i_part = 0; i_part < n_part_domains; ++i_part) {
  //   cell_center_n[i_part] = malloc( pn_cell[i_part] * sizeof(int));
  //   for(int i = 0; i < pn_cell[i_part]; ++i) {
  //     cell_center_n[i_part][i] = 1;
  //   }
  // }

  // PDM_part_to_block_exch(ptb_equi,
  //                        3 * sizeof(double),
  //                        PDM_STRIDE_VAR_INTERLACED,
  //                        1,
  //                        cell_center_n,
  //             (void **)  cell_center,
  //                        &blk_center_n,
  //              (void **) &blk_cell_center);

  // free(blk_center_n);
  // for(int i_part = 0; i_part < n_part_domains; ++i_part) {
  //   free(cell_center_n[i_part]);
  // }
  // free(cell_center_n);

  // PDM_log_trace_array_long(dextract_gnum, dn_cell_equi, "dextract_gnum : ");

  // PDM_log_trace_array_double(blk_cell_center, 3 * dn_cell_equi, "blk_cell_center : ");

  // PDM_g_num_t* distrib_extract_cell_equi = (PDM_g_num_t *) PDM_part_to_block_distrib_index_get(ptb_equi);
  // PDM_log_trace_array_long(distrib_extract_cell_equi, n_rank+1, "distrib_extract_cell_equi : ");

  // PDM_g_num_t* dequi_cell_ln_to_gn = malloc( dn_cell_equi * sizeof(PDM_g_num_t));
  // PDM_dreorder_from_coords(PDM_PART_GEOM_HILBERT,
  //                          3,
  //                          distrib_extract_cell_equi,
  //                          blk_cell_center,
  //                          dequi_cell_ln_to_gn,
  //                          comm);

  // PDM_log_trace_array_long(dequi_cell_ln_to_gn, dn_cell_equi, "dequi_cell_ln_to_gn : ");
  // free(dequi_cell_ln_to_gn);
  // free(blk_cell_center);


  PDM_part_to_part_t* ptp = PDM_part_to_part_create((const PDM_g_num_t **) pcell_ln_to_gn,
                                                    pn_cell,
                                                    n_part_domains,
                                                    (const PDM_g_num_t **) &dextract_gnum,
                                                    &dn_cell_equi,
                                                    1,
                                                    (const int         **) selected_g_num_idx,
                                                    (const PDM_g_num_t **) child_selected_g_num,
                                                    comm);

  /*
   * Extract cell_face
   */
  int         **pextract_cell_face_idx = (int         **) malloc( n_part_domains * sizeof(int         *));
  PDM_g_num_t **pextract_cell_face     = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  double      **pextract_cell_center   = (double      **) malloc( n_part_domains * sizeof(double      *));
  for (int i_part = 0; i_part < n_part_domains; i_part++){

    pextract_cell_face_idx[i_part] = malloc( (    pn_select_cell[i_part]+1) * sizeof(int   ));
    pextract_cell_center  [i_part] = malloc( (3 * pn_select_cell[i_part]  ) * sizeof(double));

    int i_extract_cell = 0;
    pextract_cell_face_idx[i_part][0] = 0;
    for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {
      for(int s_cell = selected_g_num_idx[i_part][i_cell]; s_cell < selected_g_num_idx[i_part][i_cell+1]; ++s_cell) {
        pextract_cell_face_idx[i_part][i_extract_cell+1] = pextract_cell_face_idx[i_part][i_extract_cell] + (pcell_face_idx[i_part][i_cell+1] - pcell_face_idx[i_part][i_cell]);

        pextract_cell_center[i_part][3*i_extract_cell  ] = cell_center[i_part][3*i_cell  ];
        pextract_cell_center[i_part][3*i_extract_cell+1] = cell_center[i_part][3*i_cell+1];
        pextract_cell_center[i_part][3*i_extract_cell+2] = cell_center[i_part][3*i_cell+2];
        i_extract_cell++;
      }
    }

    assert(i_extract_cell == pn_select_cell[i_part]);

    pextract_cell_face[i_part] = malloc( (pextract_cell_face_idx[i_part][pn_select_cell[i_part]]) * sizeof(PDM_g_num_t));

    int idx_write = 0;
    for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {
      for(int s_cell = selected_g_num_idx[i_part][i_cell]; s_cell < selected_g_num_idx[i_part][i_cell+1]; ++s_cell) {
        for(int idx_face = pcell_face_idx[i_part][i_cell]; idx_face < pcell_face_idx[i_part][i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS (pcell_face[i_part][idx_face])-1;
          int sgn    = PDM_SIGN(pcell_face[i_part][idx_face]);
          PDM_g_num_t g_face = pface_ln_to_gn[i_part][i_face];
          pextract_cell_face[i_part][idx_write++] = sgn * g_face;
        }
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_long(pextract_cell_face    [i_part], pextract_cell_face_idx[i_part][pn_select_cell[i_part]], "pextract_cell_face     : ");
      PDM_log_trace_array_int (pextract_cell_face_idx[i_part], pn_select_cell[i_part]+1                              , "pextract_cell_face_idx : ");
    }

  }

  /*
   *  Exchange cell_face in global numebering = cell_face + face_ln_to_gn
   */

  PDM_g_num_t* equi_extract_cell_face = NULL;

  int cst_stride = 6;
  int send_request = -1;
  int recv_request = -1;

  if (1 == 1) {


    PDM_part_to_part_issend(ptp,
                            sizeof(PDM_g_num_t),
                            cst_stride, // Hack here because HEXA
           (const void **)  pextract_cell_face,
                            100,
                            &send_request);

    equi_extract_cell_face = (PDM_g_num_t * ) malloc(cst_stride * dn_cell_equi * sizeof(PDM_g_num_t));

    PDM_part_to_part_irecv(ptp,
                           sizeof(PDM_g_num_t),
                           cst_stride, // Hack here because HEXA
                 (void **) &equi_extract_cell_face,
                           100,
                           &recv_request);

    PDM_part_to_part_issend_wait(ptp, send_request);
    PDM_part_to_part_irecv_wait (ptp, recv_request);
  }

  else {
    int  **part1_stride = malloc (sizeof(int *) * n_part_domains);
    void **part1_data = (void **) pextract_cell_face;
    int    request_exch;
    int  **part2_stride;

    for (int i = 0; i < n_part_domains; i++) {
      part1_stride[i] = malloc (sizeof(int) * pn_cell[i]);
      for (int j = 0; j < pn_cell[i]; j++) {
        part1_stride[i][j] = 6;
      }
    }

    PDM_g_num_t **_equi_extract_cell_face = NULL;

    PDM_part_to_part_iexch (ptp,
                            PDM_MPI_COMM_KIND_P2P,
                            PDM_STRIDE_VAR_INTERLACED,
                            PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2, // PDM_PART_TO_PART_DATA_DEF_ORDER_PART1
                            - 1,
                            sizeof(PDM_g_num_t),
                            (const int **) part1_stride,
                            (const void **)part1_data,
                            &part2_stride,
                            (void ***) &_equi_extract_cell_face,
                            &request_exch);

    equi_extract_cell_face = _equi_extract_cell_face[0];

    PDM_part_to_part_iexch_wait (ptp, request_exch);

    free (_equi_extract_cell_face);
    _equi_extract_cell_face = NULL;

    for (int i = 0; i < n_part_domains; i++) {
      free (part1_stride[i]);
    }
    free (part1_stride);

    free (part2_stride[0]);
    free (part2_stride);

  }

  if(0 == 1) {
    printf("cst_stride * dn_cell_equi : %d %d \n", cst_stride, dn_cell_equi);
    PDM_log_trace_array_long(equi_extract_cell_face, cst_stride * dn_cell_equi, "equi_extract_cell_face : ");
  }

  // PDM_MPI_Barrier(PDM_MPI_COMM_WORLD);
  printf("-- 1\n");
  fflush(stdout);
  // abort();

  /*
   * Exchange cell_center to post-treated
   */

  double* equi_extract_cell_center = (double * ) malloc(3 * dn_cell_equi * sizeof(double));


  PDM_part_to_part_issend(ptp,
                          3 * sizeof(double),
                          1,
               (const void **)  pextract_cell_center,
                          100,
                          &send_request);


  PDM_part_to_part_irecv(ptp,
                         3 * sizeof(double),
                         1, // Hack here because HEXA
                         (void **) &equi_extract_cell_center,
                         100,
                         &recv_request);

  PDM_part_to_part_issend_wait(ptp, send_request);
  PDM_part_to_part_irecv_wait (ptp, recv_request);

  printf("-- 2\n");
  fflush(stdout);

  /*
   * Cell center post with vtk
   */
  char filename[999];
  if (post) {
    sprintf(filename, "extract_cell_center_%3.3d.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              dn_cell_equi,
                              equi_extract_cell_center,
                              NULL, NULL);
  }

  printf("-- 3\n");
  fflush(stdout);

  /*
   * Prepare next step by descending connectivtity -> See / factorize _pconnectivity_with_local_num
   */
  int n_cell_face = cst_stride * dn_cell_equi;
  PDM_g_num_t *equi_parent_face_ln_to_gn = (PDM_g_num_t * ) malloc(n_cell_face * sizeof(PDM_g_num_t));
  int         *unique_order_face         = (int         * ) malloc(n_cell_face * sizeof(int        ));
  for(int i = 0; i < n_cell_face; ++i) {
    equi_parent_face_ln_to_gn[i] = PDM_ABS(equi_extract_cell_face[i]);
  }
  int n_extract_face = PDM_inplace_unique_long2(equi_parent_face_ln_to_gn, unique_order_face, 0, n_cell_face-1);

  int *equi_cell_face = (int *) malloc( n_cell_face * sizeof(int));
  for(int idx = 0; idx < n_cell_face; ++idx) {
    int g_sgn  = PDM_SIGN(equi_extract_cell_face[idx]);
    int l_elmt = unique_order_face[idx];
    equi_cell_face[idx] = (l_elmt + 1) * g_sgn;
  }

  if(0 == 1) {
    PDM_log_trace_array_long(equi_parent_face_ln_to_gn, n_extract_face, "equi_parent_face_ln_to_gn : ");
    PDM_log_trace_array_int (equi_cell_face           , n_cell_face   , "equi_cell_face : ");
  }

  /*
   * At this stage we have the face_ln_to_gn (parent) and we need to create the child
   */

  printf("-- 4\n");
  fflush(stdout);

  PDM_gen_gnum_t* gnum_extract_face = PDM_gnum_create(3,
                                                      1, // n_part
                                                      PDM_FALSE,
                                                      1.e-6,
                                                      comm,
                                                      PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_extract_face, 0, n_extract_face, equi_parent_face_ln_to_gn);
  PDM_gnum_compute(gnum_extract_face);

  PDM_g_num_t* extract_face_ln_to_gn = PDM_gnum_get(gnum_extract_face, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(extract_face_ln_to_gn, n_extract_face, "extract_face_ln_to_gn : ");
  }

  PDM_gnum_free(gnum_extract_face);

  int *equi_parent_face_idx = malloc( (n_extract_face + 1) * sizeof(int));
  equi_parent_face_idx[0] = 0;
  for(int i = 0; i < n_extract_face; ++i) {
    equi_parent_face_idx[i+1] = equi_parent_face_idx[i] + 1;
  }

  printf("-- 5\n");
  fflush(stdout);

  /*
   * Redo the same but with face_vtx
   *   No pb for orientation if face is not flip during partitionning process
   *   Todo it : part_to_part with :
   *       - part1 = extract partition
   *       - part2 = original partition
   *       - part1_to_part2 = equi_parent_face_ln_to_gn
   *  And exchange will be done by reverse in ptp in order to have in part1 the face_vtx connectivity
   */
  PDM_part_to_part_t* ptp_face = PDM_part_to_part_create((const PDM_g_num_t **) &extract_face_ln_to_gn,
                                                         &n_extract_face,
                                                         1,
                                                         (const PDM_g_num_t **) pface_ln_to_gn,
                                                         pn_face,
                                                         n_part_domains,
                                                         (const int         **) &equi_parent_face_idx,
                                                         (const PDM_g_num_t **) &equi_parent_face_ln_to_gn,
                                                         comm);

  int          *n_ref_face     = NULL;
  int         **ref_l_num_face = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp_face, &n_ref_face, &ref_l_num_face);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp_face, &gnum1_come_from_idx, &gnum1_come_from);

  if(0 == 1) {
    for (int i_part = 0; i_part < n_part_domains; i_part++) {
      PDM_log_trace_array_int (ref_l_num_face[i_part], n_ref_face[i_part]                                  , "ref_l_num_face      : ");
      PDM_log_trace_array_int (gnum1_come_from_idx[i_part], n_ref_face[i_part]                             , "gnum1_come_from_idx : ");
      PDM_log_trace_array_long(gnum1_come_from    [i_part], gnum1_come_from_idx[i_part][n_ref_face[i_part]], "gnum1_come_from     : ");
    }
  }

  // Creation buffer partition 2 : face_vtx + face_vtx_idx ( n_face )
  PDM_g_num_t **send_face_vtx = malloc(n_part_domains * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part_domains; ++i_part) {
    send_face_vtx[i_part] = malloc(4 * gnum1_come_from_idx[i_part][n_ref_face[i_part]] * sizeof(PDM_g_num_t));
    int idx_write = 0;
    for(int j = 0; j < n_ref_face[i_part]; ++j) {
      for(int k = gnum1_come_from_idx[i_part][j]; k < gnum1_come_from_idx[i_part][j+1]; ++k) {
        int l_face = ref_l_num_face[i_part][k]-1;
        for(int l = pface_vtx_idx[i_part][l_face]; l < pface_vtx_idx[i_part][l_face+1]; ++l) {
          send_face_vtx[i_part][idx_write++] = pvtx_ln_to_gn[i_part][pface_vtx[i_part][l]-1];
        }
      }
    }

    // printf("idx_write = %i | 4 * n_extract_face = %i \n", idx_write, 4 * n_extract_face);
    // PDM_log_trace_array_long(send_face_vtx[i_part], 4 * gnum1_come_from_idx[i_part][n_ref_face[i_part]], "send_face_vtx      : ");
  }

  printf("-- 6\n");
  fflush(stdout);

  PDM_part_to_part_reverse_issend(ptp_face,
                                  sizeof(PDM_g_num_t),
                                  4,
                 (const void **)  send_face_vtx, // Donc in same order of part2_to_part1_data
                                  100,
                                  &send_request);

  PDM_g_num_t *equi_extract_face_vtx = malloc( 4 * n_extract_face * sizeof(PDM_g_num_t));

  for(int i = 0; i < 4 * n_extract_face; ++i) {
    equi_extract_face_vtx[i] = -1;
  }

  PDM_part_to_part_reverse_irecv(ptp_face,
                                 sizeof(PDM_g_num_t),
                                 4, // Hack here because HEXA
                      (void **)  &equi_extract_face_vtx,  // order given by gnum1_come_from and ref_lnum2 arrays
                                 100,
                                 &recv_request);

  PDM_part_to_part_reverse_issend_wait(ptp_face, send_request);
  PDM_part_to_part_reverse_irecv_wait (ptp_face, recv_request);

  printf("-- 7\n");
  fflush(stdout);

  // PDM_log_trace_array_long(equi_extract_face_vtx, 4 * n_extract_face, "equi_extract_face_vtx : ");

  // DEBUG PART_TO_PART ONLY
  // PDM_g_num_t **send_face_vtx = malloc(n_part_domains * sizeof(PDM_g_num_t *));
  // for(int i_part = 0; i_part < n_part_domains; ++i_part) {
  //   send_face_vtx[i_part] = malloc(gnum1_come_from_idx[i_part][n_ref_face[i_part]] * sizeof(PDM_g_num_t));


  //   int idx_write = 0;
  //   for(int j = 0; j < n_ref_face[i_part]; ++j) {
  //     for(int k = gnum1_come_from_idx[i_part][j]; k < gnum1_come_from_idx[i_part][j+1]; ++k) {
  //       int l_face = ref_l_num_face[i_part][k]-1;
  //       send_face_vtx[i_part][idx_write++] = pface_ln_to_gn[i_part][l_face];
  //     }
  //   }

  //   printf("idx_write = %i | 4 * n_extract_face = %i \n", idx_write, 4 * n_extract_face);
  //   PDM_log_trace_array_long(send_face_vtx[i_part], 1 * gnum1_come_from_idx[i_part][n_ref_face[i_part]], "send_face_vtx      : ");
  // }


  // PDM_part_to_part_reverse_issend(ptp_face,
  //                                 sizeof(PDM_g_num_t),
  //                                 1,
  //                      (void **)  send_face_vtx, // Donc in same order of part2_to_part1_data
  //                                 100,
  //                                 &send_request);

  // PDM_g_num_t *equi_extract_face_vtx = malloc( 1 * n_extract_face * sizeof(PDM_g_num_t));

  // for(int i = 0; i < 1 * n_extract_face; ++i) {
  //   equi_extract_face_vtx[i] = -1;
  // }

  // PDM_part_to_part_reverse_irecv(ptp_face,
  //                                sizeof(PDM_g_num_t),
  //                                1, // Hack here because HEXA
  //                     (void **)  &equi_extract_face_vtx,  // order given by gnum1_come_from and ref_lnum2 arrays
  //                                100,
  //                                &recv_request);

  // PDM_part_to_part_reverse_issend_wait(ptp_face, send_request);
  // PDM_part_to_part_reverse_irecv_wait (ptp_face, recv_request);
  // PDM_log_trace_array_long(equi_extract_face_vtx, 1 * n_extract_face, "equi_extract_face_vtx : ");

  PDM_part_to_part_free(ptp_face);

  /*
   * Redo the same for vtx
   */
  int n_face_vtx = 4 * n_extract_face;
  PDM_g_num_t *equi_parent_vtx_ln_to_gn = (PDM_g_num_t * ) malloc(n_face_vtx * sizeof(PDM_g_num_t));
  int         *unique_order_vtx         = (int         * ) malloc(n_face_vtx * sizeof(int        ));
  for(int i = 0; i < n_face_vtx; ++i) {
    equi_parent_vtx_ln_to_gn[i] = PDM_ABS(equi_extract_face_vtx[i]);
  }
  int n_extract_vtx = PDM_inplace_unique_long2(equi_parent_vtx_ln_to_gn, unique_order_vtx, 0, n_face_vtx-1);

  int *equi_face_vtx = (int *) malloc( n_face_vtx * sizeof(int));
  for(int idx = 0; idx < n_face_vtx; ++idx) {
    int g_sgn  = PDM_SIGN(equi_extract_face_vtx[idx]);
    int l_elmt = unique_order_vtx[idx];
    equi_face_vtx[idx] = (l_elmt + 1) * g_sgn;
  }
  free(unique_order_vtx);

  if(0 == 1) {
    PDM_log_trace_array_long(equi_parent_vtx_ln_to_gn, n_extract_vtx, "equi_parent_vtx_ln_to_gn : ");
    PDM_log_trace_array_int (equi_face_vtx           , n_face_vtx   , "equi_face_vtx            : ");
  }

  // TO see with Eric, tag pb because cast into int
  int          *part2_n_extract_vtx        = NULL;
  int         **part2_face_vtx_idx         = NULL;
  int         **part2_face_vtx             = NULL;
  PDM_g_num_t **part2_parent_vtx_ln_to_gn  = NULL;
  PDM_g_num_t **part2_vtx_ln_to_gn         = NULL;
  PDM_pconnectivity_to_pconnectivity(comm,
                                     n_part,
            (const int            *) pn_face,
            (const int           **) pface_vtx_idx,
            (const int           **) pface_vtx,
            (const PDM_g_num_t   **) pface_ln_to_gn,
            (const PDM_g_num_t   **) pvtx_ln_to_gn,
                                     1,
            (const int            *) &n_extract_face,
            (const PDM_g_num_t   **) &extract_face_ln_to_gn,
            (const int           **) &equi_parent_face_idx,
            (const PDM_g_num_t   **) &equi_parent_face_ln_to_gn,
                                     &part2_n_extract_vtx,
                                     &part2_face_vtx_idx,
                                     &part2_face_vtx,
                                     &part2_parent_vtx_ln_to_gn,
                                     &part2_vtx_ln_to_gn);
  free(part2_n_extract_vtx);
  free(part2_face_vtx_idx[0]);
  free(part2_face_vtx[0]);
  free(part2_vtx_ln_to_gn[0]);
  free(part2_parent_vtx_ln_to_gn[0]);
  free(part2_face_vtx_idx);
  free(part2_face_vtx);
  free(part2_vtx_ln_to_gn);
  free(part2_parent_vtx_ln_to_gn);



  /*
   * At this stage we have the vtx_ln_to_gn (parent) and we need to create the child
   */
  PDM_gen_gnum_t* gnum_extract_vtx = PDM_gnum_create(3,
                                                      1, // n_part
                                                      PDM_FALSE,
                                                      1.e-6,
                                                      comm,
                                                      PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_extract_vtx, 0, n_extract_vtx, equi_parent_vtx_ln_to_gn);
  PDM_gnum_compute(gnum_extract_vtx);

  PDM_g_num_t* extract_vtx_ln_to_gn = PDM_gnum_get(gnum_extract_vtx, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(extract_vtx_ln_to_gn, n_extract_vtx, "extract_vtx_ln_to_gn : ");
  }

  PDM_gnum_free(gnum_extract_vtx);

  int *equi_parent_vtx_idx = malloc( (n_extract_vtx + 1) * sizeof(int));
  equi_parent_vtx_idx[0] = 0;
  for(int i = 0; i < n_extract_vtx; ++i) {
    equi_parent_vtx_idx[i+1] = equi_parent_vtx_idx[i] + 1;
  }

  printf("-- 8\n");
  fflush(stdout);

  PDM_part_to_part_t* ptp_vtx = PDM_part_to_part_create((const PDM_g_num_t **) &extract_vtx_ln_to_gn,
                                                         &n_extract_vtx,
                                                         1,
                                                         (const PDM_g_num_t **) pvtx_ln_to_gn,
                                                         pn_vtx,
                                                         n_part_domains,
                                                         (const int         **) &equi_parent_vtx_idx,
                                                         (const PDM_g_num_t **) &equi_parent_vtx_ln_to_gn,
                                                         comm);

  int          *n_ref_vtx     = NULL;
  int         **ref_l_num_vtx = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp_vtx, &n_ref_vtx, &ref_l_num_vtx);

  int         **gnum1_come_from_vtx_idx = NULL;
  PDM_g_num_t **gnum1_come_from_vtx     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp_vtx, &gnum1_come_from_vtx_idx, &gnum1_come_from_vtx);

  if(0 == 1) {
    for (int i_part = 0; i_part < n_part_domains; i_part++) {
      PDM_log_trace_array_int (ref_l_num_vtx      [i_part], n_ref_vtx[i_part]                                     , "ref_l_num_vtx           : ");
      PDM_log_trace_array_int (gnum1_come_from_vtx_idx[i_part], n_ref_vtx[i_part]                                 , "gnum1_come_from_vtx_idx : ");
      PDM_log_trace_array_long(gnum1_come_from_vtx    [i_part], gnum1_come_from_vtx_idx[i_part][n_ref_vtx[i_part]], "gnum1_come_from_vtx     : ");
    }
  }

  // Creation buffer partition 2 : vtx_coord ( n_face )
  double **send_vtx_coord = malloc(n_part_domains * sizeof(double *));
  for(int i_part = 0; i_part < n_part_domains; ++i_part) {
    send_vtx_coord[i_part] = malloc(3 * gnum1_come_from_vtx_idx[i_part][n_ref_vtx[i_part]] * sizeof(double));
    int idx_write = 0;
    for(int j = 0; j < n_ref_vtx[i_part]; ++j) {
      for(int k = gnum1_come_from_vtx_idx[i_part][j]; k < gnum1_come_from_vtx_idx[i_part][j+1]; ++k) {
        int l_vtx = ref_l_num_vtx[i_part][k]-1;
        send_vtx_coord[i_part][3*idx_write  ] = pvtx_coord[i_part][3*l_vtx  ];
        send_vtx_coord[i_part][3*idx_write+1] = pvtx_coord[i_part][3*l_vtx+1];
        send_vtx_coord[i_part][3*idx_write+2] = pvtx_coord[i_part][3*l_vtx+2];
        idx_write++;
      }
    }

    // PDM_log_trace_array_double(send_vtx_coord[i_part], 3 * gnum1_come_from_vtx_idx[i_part][n_ref_vtx[i_part]], "send_vtx_coord      : ");
  }

  double *equi_extract_vtx_coord = NULL;

  if (1 == 1) {

    PDM_part_to_part_reverse_issend(ptp_vtx,
                                    sizeof(double),
                                    3,
                   (const void **)  send_vtx_coord, // Donc in same order of part2_to_part1_data
                                    100,
                                    &send_request);


    equi_extract_vtx_coord = malloc( 3 * n_extract_vtx * sizeof(double));

    for(int i = 0; i < 3 * n_extract_vtx; ++i) {
      equi_extract_vtx_coord[i] = -1;
    }

    PDM_part_to_part_reverse_irecv(ptp_vtx,
                                   sizeof(double),
                                   3,
                        (void **)  &equi_extract_vtx_coord,  // order given by gnum1_come_from and ref_lnum2 arrays
                                   100,
                                   &recv_request);

    PDM_part_to_part_reverse_issend_wait(ptp_vtx, send_request);
    PDM_part_to_part_reverse_irecv_wait (ptp_vtx, recv_request);

  }

  else {

    int      request_exch;

    int    **part2_stride = malloc (sizeof(int*) * n_part_domains);

    for(int i_part = 0; i_part < n_part_domains; ++i_part) {
      part2_stride[i_part] = malloc(gnum1_come_from_vtx_idx[i_part][n_ref_vtx[i_part]] * sizeof(int));
      for (int j = 0; j < gnum1_come_from_vtx_idx[i_part][n_ref_vtx[i_part]]; j++) {
        part2_stride[i_part][j] = 3;
      }
    }

    int    **part1_stride = NULL;;
    double **part1_data;

    PDM_part_to_part_reverse_iexch (ptp_vtx,
                                    PDM_MPI_COMM_KIND_P2P,
                                    PDM_STRIDE_VAR_INTERLACED,
                                    PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                    -1,
                                    sizeof(double),
                   (const int **)   part2_stride,
                   (const void **)  send_vtx_coord,
                                   &part1_stride,
                   (void ***)      &part1_data,
                                   &request_exch);

    PDM_part_to_part_reverse_iexch_wait (ptp_vtx, request_exch);

    equi_extract_vtx_coord = part1_data[0];
    free (part1_data);
    free (part1_stride[0]);
    free (part1_stride);
    for(int i_part = 0; i_part < n_part_domains; ++i_part) {
      free (part2_stride[i_part]);
    }
    free (part2_stride);
  }

  PDM_part_to_part_free(ptp_vtx);

  printf("-- 9\n");
  fflush(stdout);



  if (post) {
    PDM_log_trace_array_double(equi_extract_vtx_coord, 3 * n_extract_vtx, "equi_extract_vtx_coord : ");

    sprintf(filename, "extract_vtx_coord_%3.3d.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_extract_vtx,
                              equi_extract_vtx_coord,
                              NULL, NULL);
  }

  int* equi_face_vtx_idx = malloc( (n_extract_face + 1) * sizeof(int));
  equi_face_vtx_idx[0] = 0;
  for(int i = 0; i < n_extract_face; ++i) {
    equi_face_vtx_idx[i+1] = equi_face_vtx_idx[i] + 4;
  }

  if (post) {
    sprintf(filename, "face_vtx_coord_%3.3d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           n_extract_vtx,
                           equi_extract_vtx_coord,
                           extract_vtx_ln_to_gn,
                           n_extract_face,
                           equi_face_vtx_idx,
                           equi_face_vtx,
                           extract_face_ln_to_gn,
                           NULL);
  }

  free(equi_face_vtx_idx);


  free(equi_extract_vtx_coord);

  free(equi_parent_face_ln_to_gn);
  free(equi_parent_vtx_ln_to_gn);
  free(extract_face_ln_to_gn);
  free(unique_order_face);
  free(equi_cell_face);
  free(equi_face_vtx);
  free(equi_parent_face_idx);
  free(equi_extract_face_vtx);
  free(equi_parent_vtx_idx);


  free(equi_extract_cell_face);
  free(equi_extract_cell_center);
  free(extract_vtx_ln_to_gn);

  PDM_part_to_part_free(ptp);

  PDM_part_to_block_free(ptb_equi);
  PDM_gnum_free(gnum_extract);


  for (int i_part = 0; i_part < n_part_domains; i_part++){
    free(cell_center       [i_part]);
    free(selected_g_num    [i_part]);
    free(selected_g_num_idx[i_part]);
    free(weight[i_part]);
    free(tmp_extract_cell_center[i_part]);
    free(pextract_cell_face[i_part]);
    free(pextract_cell_face_idx[i_part]);
    free(pextract_cell_center[i_part]);
    free(send_face_vtx[i_part]);
    free(send_vtx_coord[i_part]);
  }
  free(cell_center);
  free(selected_g_num);
  free(selected_g_num_idx);
  free(pn_cell);
  free(pn_face);
  free(pn_vtx);
  free(send_face_vtx);
  free(send_vtx_coord);
  free(pcell_face);
  free(pcell_face_idx);
  free(pvtx_coord);
  free(pface_vtx);
  free(pface_vtx_idx);
  free(child_selected_g_num);
  free(pn_select_cell);
  free(weight);
  free(tmp_extract_cell_center);
  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pextract_cell_face);
  free(pextract_cell_face_idx);
  free(pextract_cell_center);

  PDM_multipart_free(mpart);
  PDM_dcube_gen_free(dcube);
  PDM_dmesh_free(dm);


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize();

  return 0;
}
