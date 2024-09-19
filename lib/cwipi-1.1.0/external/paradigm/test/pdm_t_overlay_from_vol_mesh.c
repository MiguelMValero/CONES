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
#include "pdm_extract_part.h"
#include "pdm_vtk.h"
#include "pdm_dmesh.h"
#include "pdm_unique.h"
#include "pdm_part_geom.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_overlay.h"

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
           PDM_g_num_t  *n_surf_vtx_seg,
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
    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_surf_vtx_seg = atol(argv[i]);
        *n_surf_vtx_seg = (PDM_g_num_t) _n_surf_vtx_seg;
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

static
void
_generate_surf_mesh_1d
(
 PDM_MPI_Comm    comm,
 double          length,
 double          zero_x,
 double          zero_y,
 double          zero_z,
 PDM_g_num_t     surf_nx_vtx_seg,
 PDM_g_num_t     surf_ny_vtx_seg,
 int            *psurf_n_face,
 int            *psurf_n_vtx,
 int           **surf_face_vtx_idx,
 int           **surf_face_vtx,
 PDM_g_num_t   **surf_face_ln_to_gn,
 PDM_g_num_t   **surf_vtx_ln_to_gn,
 double        **surf_vtx_coord
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t gn_vtx  = (surf_nx_vtx_seg    ) * (surf_ny_vtx_seg    );
  PDM_g_num_t gn_face = (surf_nx_vtx_seg - 1) * (surf_ny_vtx_seg - 1);

  int dcube_nx = surf_nx_vtx_seg - 1;

  PDM_g_num_t* distrib_face = PDM_compute_uniform_entity_distribution(comm, gn_face);
  PDM_g_num_t* distrib_vtx  = PDM_compute_uniform_entity_distribution(comm, gn_vtx);

  int dn_vtx  = (int) (distrib_vtx[i_rank+1] - distrib_vtx[i_rank]);
  int dn_face = (int) (distrib_face[i_rank+1] - distrib_face[i_rank]);

  double *dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);

  double step_x = length / (double) (surf_nx_vtx_seg - 1);
  double step_y = length / (double) (surf_ny_vtx_seg - 1);

  for (int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {

    PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

    PDM_g_num_t indi = g_vtx % surf_nx_vtx_seg;
    PDM_g_num_t indj = g_vtx / surf_nx_vtx_seg;

    dvtx_coord[3 * i_vtx    ] = zero_x;
    dvtx_coord[3 * i_vtx + 1] = indi * step_x + zero_y;
    dvtx_coord[3 * i_vtx + 2] = indj * step_y + zero_z;
  }

  int         *dface_vtx_idx = malloc( (dn_face+1) * sizeof(int        ));
  PDM_g_num_t *dface_vtx     = malloc( 4 * dn_face * sizeof(PDM_g_num_t));

  dface_vtx_idx[0] = 0;
  for (int i_face = 0; i_face < dn_face; ++i_face) {

    dface_vtx_idx[i_face+1] = dface_vtx_idx[i_face] + 4;
    PDM_g_num_t g = distrib_face[i_rank] + i_face;

    PDM_g_num_t indi = g % dcube_nx;
    PDM_g_num_t indj = g / dcube_nx;

    dface_vtx[dface_vtx_idx[i_face]  ] = 1 + (indi  ) + (dcube_nx + 1) * (indj    );
    dface_vtx[dface_vtx_idx[i_face]+1] = 1 + (indi+1) + (dcube_nx + 1) * (indj    );
    dface_vtx[dface_vtx_idx[i_face]+2] = 1 + (indi+1) + (dcube_nx + 1) * (indj + 1);
    dface_vtx[dface_vtx_idx[i_face]+3] = 1 + (indi  ) + (dcube_nx + 1) * (indj + 1);


  }

  // PDM_log_trace_connectivity_long(dface_vtx_idx, dface_vtx, dn_face, "dface_vtx");

  int          _psurf_n_face       = dn_face;
  PDM_g_num_t *_surf_face_ln_to_gn = malloc(dn_face * sizeof(PDM_g_num_t));
  for(int i_face = 0; i_face < dn_face; ++i_face) {
    _surf_face_ln_to_gn[i_face] = distrib_face[i_rank] + i_face + 1;
  }

  int          _psurf_n_vtx        = 0;
  PDM_g_num_t *_surf_vtx_ln_to_gn  = NULL;
  int         *_surf_face_vtx_idx  = NULL;
  int         *_surf_face_vtx      = NULL;
  double      *_surf_vtx_coord     = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           _psurf_n_face,
                                                           _surf_face_ln_to_gn,
                                                           &_psurf_n_vtx,
                                                           &_surf_vtx_ln_to_gn,
                                                           &_surf_face_vtx_idx,
                                                           &_surf_face_vtx);

  double **tmp_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_vtx,
                                        dvtx_coord,
                                        &_psurf_n_vtx,
                 (const PDM_g_num_t **) &_surf_vtx_ln_to_gn,
                                        &tmp_vtx_coord);
  _surf_vtx_coord = tmp_vtx_coord[0];
  free(tmp_vtx_coord);


  *psurf_n_face       = _psurf_n_face;
  *psurf_n_vtx        = _psurf_n_vtx;
  *surf_face_ln_to_gn = _surf_face_ln_to_gn;
  *surf_vtx_ln_to_gn  = _surf_vtx_ln_to_gn;
  *surf_face_vtx_idx  = _surf_face_vtx_idx;
  *surf_face_vtx      = _surf_face_vtx;
  *surf_vtx_coord     = _surf_vtx_coord;


  free(dvtx_coord);
  free(distrib_face);
  free(distrib_vtx);
  free(dface_vtx_idx);
  free(dface_vtx);
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

  PDM_g_num_t        n_vtx_seg      = 10;
  PDM_g_num_t        n_surf_vtx_seg = 10;
  double             length         = 1.;
  int                n_part         = 1;
  int                post           = 0;

  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &n_surf_vtx_seg,
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

  int         **selected_face_l_num     = (int         **) malloc( n_part_domains * sizeof(int         *));
  PDM_g_num_t **pcell_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pface_ln_to_gn          = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pvtx_ln_to_gn           = (PDM_g_num_t **) malloc( n_part_domains * sizeof(PDM_g_num_t *));
  int          *pn_cell                 = (int          *) malloc( n_part_domains * sizeof(int          ));
  int          *pn_face                 = (int          *) malloc( n_part_domains * sizeof(int          ));
  int          *pn_vtx                  = (int          *) malloc( n_part_domains * sizeof(int          ));
  int          *pn_select_face          = (int          *) malloc( n_part_domains * sizeof(int          ));
  // double      **weight                  = (double      **) malloc( n_part_domains * sizeof(double      *));
  int         **pcell_face              = (int         **) malloc( n_part_domains * sizeof(int         *));
  int         **pcell_face_idx          = (int         **) malloc( n_part_domains * sizeof(int         *));
  int         **pface_vtx               = (int         **) malloc( n_part_domains * sizeof(int         *));
  int         **pface_vtx_idx           = (int         **) malloc( n_part_domains * sizeof(int         *));
  double      **pvtx_coord              = (double      **) malloc( n_part_domains * sizeof(double      *));

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

    int          n_group           = 0;
    int*         group_face_idx    = NULL;
    int*         group_face        = NULL;
    PDM_g_num_t* group_face_lngn   = NULL;
    PDM_multipart_group_get(mpart,
                            0,
                            i_part,
                            PDM_MESH_ENTITY_FACE,
                            &n_group,
                            &group_face_idx,
                            &group_face,
                            &group_face_lngn,
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

    selected_face_l_num         [i_part] = (int         *) malloc(  n_cell          * sizeof(int        ));

    /*
     * Sub-part
     */
    int i_extract_face_group = 2;
    int n_select_face        = 0;

    for(int idx_face = group_face_idx[i_extract_face_group]; idx_face < group_face_idx[i_extract_face_group+1]; ++idx_face) {
      int i_face = group_face[idx_face] - 1;
      selected_face_l_num[i_part][n_select_face] = i_face;
      n_select_face++;
    }

    selected_face_l_num[i_part] = realloc(selected_face_l_num[i_part], n_select_face * sizeof(int        ));
    pn_select_face     [i_part] = n_select_face;

    // PDM_log_trace_array_int (selected_face_l_num[i_part], n_select_face, "selected_face_l_num : ");
    free(face_center);

  }

  /*
   * Extract
   */
  int n_part_out = 1;
  PDM_extract_part_t* extrp = PDM_extract_part_create(2,
                                                      n_part,
                                                      n_part_out,
                                                      PDM_EXTRACT_PART_KIND_LOCAL,
                                                      PDM_SPLIT_DUAL_WITH_HILBERT,
                                                      PDM_TRUE, // compute_child_gnum
                                                      PDM_OWNERSHIP_KEEP,
                                                      comm);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_extract_part_part_set(extrp,
                              i_part,
                              pn_cell[i_part],
                              pn_face[i_part],
                              -1, // pn_edge[i_part],
                              pn_vtx[i_part],
                              pcell_face_idx[i_part],
                              pcell_face[i_part],
                              NULL, //pface_edge_idx[i_part],
                              NULL, //pface_edge[i_part],
                              NULL, //pedge_vtx[i_part],
                              pface_vtx_idx[i_part],
                              pface_vtx[i_part],
                              pcell_ln_to_gn[i_part],
                              pface_ln_to_gn[i_part],
                              NULL, //pedge_ln_to_gn[i_part],
                              pvtx_ln_to_gn[i_part],
                              pvtx_coord[i_part]);

    PDM_extract_part_selected_lnum_set(extrp,
                                       i_part,
                                       pn_select_face[i_part],
                                       selected_face_l_num[i_part]);

    // PDM_log_trace_array_int(selected_face_l_num[i_part], pn_select_face[i_part], "selected_face_l_num ::");

  }


  PDM_extract_part_compute(extrp);

  int          *pn_extract_face        = malloc(n_part_out * sizeof(int          ));
  int          *pn_extract_vtx         = malloc(n_part_out * sizeof(int          ));
  int         **pextract_face_vtx      = malloc(n_part_out * sizeof(int         *));
  int         **pextract_face_vtx_idx  = malloc(n_part_out * sizeof(int         *));
  double      **pextract_vtx           = malloc(n_part_out * sizeof(double      *));
  PDM_g_num_t **pextract_face_ln_to_gn = malloc(n_part_out * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pextract_vtx_ln_to_gn  = malloc(n_part_out * sizeof(PDM_g_num_t *));


  for(int i_part = 0; i_part < n_part_out; ++i_part) {

    pn_extract_face[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE);

    pn_extract_vtx[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                           i_part,
                                                           PDM_MESH_ENTITY_VTX);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                      &pextract_face_vtx[i_part],
                                      &pextract_face_vtx_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    PDM_extract_part_vtx_coord_get(extrp,
                                   i_part,
                                   &pextract_vtx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &pextract_face_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_VTX,
                                  &pextract_vtx_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);
  }

  /*
   *  Create simple mesh
   */
  int           psurf_n_face       = 0;
  int           psurf_n_vtx        = 0;
  int          *surf_face_vtx_idx  = NULL;
  int          *surf_face_vtx      = NULL;
  PDM_g_num_t  *surf_face_ln_to_gn = NULL;
  PDM_g_num_t  *surf_vtx_ln_to_gn  = NULL;
  double       *surf_vtx_coord     = NULL;
  _generate_surf_mesh_1d(comm,
                         length,
                         0.,
                         0.,
                         0.,
                         n_surf_vtx_seg,
                         2,
                         &psurf_n_face,
                         &psurf_n_vtx,
                         &surf_face_vtx_idx,
                         &surf_face_vtx,
                         &surf_face_ln_to_gn,
                         &surf_vtx_ln_to_gn,
                         &surf_vtx_coord);

  /*
   * Export vtk en légende
   */
  if(post) {
    char filename[999];
    for(int i_part = 0; i_part < n_part_out; ++i_part) {

      sprintf(filename, "extract_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                pn_extract_vtx[i_part],
                                pextract_vtx[i_part],
                                NULL, NULL);

      sprintf(filename, "extract_face_vtx_coord_%3.3d_%3.3d.vtk", i_part, i_rank);
      PDM_vtk_write_polydata(filename,
                             pn_extract_vtx[i_part],
                             pextract_vtx[i_part],
                             pextract_vtx_ln_to_gn[i_part],
                             pn_extract_face[i_part],
                             pextract_face_vtx_idx[i_part],
                             pextract_face_vtx[i_part],
                             pextract_face_ln_to_gn[i_part],
                             NULL);

    }

    sprintf(filename, "extract_surf_face_vtx_coord_%3.3d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           psurf_n_vtx,
                           surf_vtx_coord,
                           surf_vtx_ln_to_gn,
                           psurf_n_face,
                           surf_face_vtx_idx,
                           surf_face_vtx,
                           surf_face_ln_to_gn,
                           NULL);
  }


  /*
   * Prepare overlay
   */
  double projectCoeff = 0.;
  int n_part2 = 1;
  PDM_ol_t *ol = PDM_ol_create (n_part,
                                n_part2,
                                projectCoeff,
                                PDM_MPI_COMM_WORLD);

  PDM_ol_parameter_set (ol,
                        PDM_OL_CAR_LENGTH_TOL,
                        1e-4);

  PDM_ol_parameter_set (ol,
                        PDM_OL_EXTENTS_TOL,
                        1e-4);

  /* Set mesh A */
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_ol_input_mesh_set (ol,
                           PDM_OL_MESH_A,
                           i_part,
                           pn_extract_face       [i_part],
                           pextract_face_vtx_idx [i_part],
                           pextract_face_vtx     [i_part],
                           pextract_face_ln_to_gn[i_part],
                           pn_extract_vtx        [i_part],
                           pextract_vtx          [i_part],
                           pextract_vtx_ln_to_gn [i_part]);
  }

  /* Set mesh B */
  PDM_ol_input_mesh_set (ol,
                         PDM_OL_MESH_B,
                         0,
                         psurf_n_face,
                         surf_face_vtx_idx,
                         surf_face_vtx,
                         surf_face_ln_to_gn,
                         psurf_n_vtx,
                         surf_vtx_coord,
                         surf_vtx_ln_to_gn);

  PDM_ol_compute (ol);

  if (i_rank == 0){
    PDM_ol_dump_times (ol);
  }

  PDM_ol_mesh_t mesht = PDM_OL_MESH_A;
  for(int i_mesh = 0; i_mesh < 2; ++i_mesh) {
    for (int ipart = 0; ipart < n_part; ipart++) {

      int           n_ol_face;
      int           n_ol_linked_face;
      int           n_ol_vtx;
      int           s_ol_face_ini_vtx;
      int           s_ol_face_vtx;
      int           s_init_to_ol_face;

      PDM_ol_part_mesh_dim_get (ol,
                                mesht,
                                ipart,
                                &n_ol_face,
                                &n_ol_linked_face,
                                &n_ol_vtx,
                                &s_ol_face_ini_vtx,
                                &s_ol_face_vtx,
                                &s_init_to_ol_face);


      int            *ol_face_ini_vtx_idx;
      int            *ol_face_ini_vtx;
      int            *ol_face_vtx_idx;
      int            *ol_face_vtx;
      int            *ol_linked_face_proc_idx;
      int            *ol_linked_face;
      PDM_g_num_t    *ol_face_ln_to_gn;
      double         *ol_vtx_coords;
      PDM_g_num_t    *ol_vtx_ln_to_gn;
      int            *init_to_ol_face_idx;
      int            *init_to_ol_face;

      PDM_ol_mesh_entities_get (ol,
                                mesht,
                                ipart,
                                &ol_face_ini_vtx_idx,
                                &ol_face_ini_vtx,
                                &ol_face_vtx_idx,
                                &ol_face_vtx,
                                &ol_linked_face_proc_idx,
                                &ol_linked_face,
                                &ol_face_ln_to_gn,
                                &ol_vtx_coords,
                                &ol_vtx_ln_to_gn,
                                &init_to_ol_face_idx,
                                &init_to_ol_face);

      if (post) {
        char filename[999];
        if(mesht == PDM_OL_MESH_A) {
          sprintf(filename, "mesh_a_inter_b_%3.3d.vtk", i_rank);
        } else {
          sprintf(filename, "mesh_b_inter_a_%3.3d.vtk", i_rank);
        }
        PDM_vtk_write_polydata(filename,
                               n_ol_vtx,
                               ol_vtx_coords,
                               ol_vtx_ln_to_gn,
                               n_ol_face,
                               ol_face_vtx_idx,
                               ol_face_vtx,
                               ol_face_ln_to_gn,
                               NULL);
      }

      free(ol_face_ini_vtx_idx);
      free(ol_face_ini_vtx);
      free(ol_face_vtx_idx);
      free(ol_face_vtx);
      free(ol_linked_face_proc_idx);
      free(ol_linked_face);
      free(ol_face_ln_to_gn);
      free(ol_vtx_coords);
      free(ol_vtx_ln_to_gn);
      free(init_to_ol_face_idx);
      free(init_to_ol_face);

    }
    mesht = PDM_OL_MESH_B;
  }




  PDM_ol_del (ol);
  PDM_extract_part_free(extrp);

  free(surf_face_vtx_idx  );
  free(surf_face_vtx      );
  free(surf_face_ln_to_gn );
  free(surf_vtx_ln_to_gn  );
  free(surf_vtx_coord     );

  for (int i_part = 0; i_part < n_part_domains; i_part++){
    free(selected_face_l_num[i_part]);
  }
  free(selected_face_l_num);
  free(pn_cell);
  free(pn_face);
  free(pn_vtx);
  free(pn_select_face);

  free(pcell_ln_to_gn);
  free(pface_ln_to_gn);
  free(pvtx_ln_to_gn );
  free(pcell_face    );
  free(pcell_face_idx);
  free(pface_vtx     );
  free(pface_vtx_idx );
  free(pvtx_coord    );

  free(pn_extract_face       );
  free(pn_extract_vtx        );
  free(pextract_face_vtx     );
  free(pextract_face_vtx_idx );
  free(pextract_vtx          );
  free(pextract_face_ln_to_gn);
  free(pextract_vtx_ln_to_gn );

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
