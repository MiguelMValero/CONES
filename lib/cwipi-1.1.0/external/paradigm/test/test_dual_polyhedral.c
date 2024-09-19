#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_para_graph_dual.h"
#include "pdm_part.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_dmesh.h"
#include "pdm_printf.h"
#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_geom_elem.h"

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
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int           argc,
           char        **argv,
           PDM_g_num_t  *nx,
           PDM_g_num_t  *ny,
           PDM_g_num_t  *nz,
           double       *lengthx,
           double       *lengthy,
           double       *lengthz,
           int          *n_part,
           int          *randomize,
           int          *random_seed,
           int          *post,
           int          *method)
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
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
        *ny = (PDM_g_num_t) _n;
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *ny = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthx = atof(argv[i]);
        *lengthy = atof(argv[i]);
        *lengthz = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-lx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthx = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ly") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthy = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-lz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthz = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      *randomize = 1;
    }
    else if (strcmp(argv[i], "-seed") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *random_seed = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static void
_compute_dual_closed
(
 int          n_cell,
 int          n_face,
 int          n_edge,
 int          n_vtx,
 int         *cell_face_idx,
 PDM_g_num_t *cell_face,
 int         *face_edge_idx,
 PDM_g_num_t *face_edge,
 int         *edge_vtx_idx,
 PDM_g_num_t *edge_vtx,
 double      *coords
)
{
  double* edge_center = (double *) malloc(3 * n_edge * sizeof(double));
  double* face_center = (double *) malloc(3 * n_face * sizeof(double));
  double* cell_center = (double *) malloc(3 * n_cell * sizeof(double));
  double* edge_surf   = (double *) malloc(3 * n_edge * sizeof(double));

  for(int i_edge = 0; i_edge < n_edge; ++i_edge) {
    int i_vtx1 = edge_vtx[2*i_edge  ]-1;
    int i_vtx2 = edge_vtx[2*i_edge+1]-1;
    edge_center[3*i_edge  ] = 0.5 * (coords[3*i_vtx1  ] + coords[3*i_vtx2  ]);
    edge_center[3*i_edge+1] = 0.5 * (coords[3*i_vtx1+1] + coords[3*i_vtx2+1]);
    edge_center[3*i_edge+2] = 0.5 * (coords[3*i_vtx1+2] + coords[3*i_vtx2+2]);
    edge_surf[3*i_edge  ] = 0.;
    edge_surf[3*i_edge+1] = 0.;
    edge_surf[3*i_edge+2] = 0.;

    // printf("x1, y1, z1 = %12.5e %12.5e %12.5e | x2, y2, z2 = %12.5e %12.5e %12.5e \n", coords[3*i_vtx1  ],
    //                                                                                    coords[3*i_vtx1+1],
    //                                                                                    coords[3*i_vtx1+2],
    //                                                                                    coords[3*i_vtx2  ],
    //                                                                                    coords[3*i_vtx2+1],
    //                                                                                    coords[3*i_vtx2+2]);
    // printf("xe, ye, ze = %12.5e %12.5e %12.5e \n", edge_center[3*i_edge  ], edge_center[3*i_edge+1], edge_center[3*i_edge+2]);
  }
  // exit(1);

  for(int i_face = 0; i_face < n_face; ++i_face) {

    face_center[3*i_face  ] = 0.;
    face_center[3*i_face+1] = 0.;
    face_center[3*i_face+2] = 0.;
    double pond = 1./(face_edge_idx[i_face+1] - face_edge_idx[i_face]);

    for(int idx_edge = face_edge_idx[i_face]; idx_edge < face_edge_idx[i_face+1]; ++idx_edge) {
      int i_edge = PDM_ABS(face_edge[idx_edge])-1;
      face_center[3*i_face  ] += edge_center[3*i_edge  ];
      face_center[3*i_face+1] += edge_center[3*i_edge+1];
      face_center[3*i_face+2] += edge_center[3*i_edge+2];
    }
    face_center[3*i_face  ] = face_center[3*i_face  ] * pond;
    face_center[3*i_face+1] = face_center[3*i_face+1] * pond;
    face_center[3*i_face+2] = face_center[3*i_face+2] * pond;
  }

  for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

    cell_center[3*i_cell  ] = 0.;
    cell_center[3*i_cell+1] = 0.;
    cell_center[3*i_cell+2] = 0.;
    double pond = 1./(cell_face_idx[i_cell+1] - cell_face_idx[i_cell]);

    for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {
      int i_face = PDM_ABS(cell_face[idx_face])-1;
      cell_center[3*i_cell  ] += face_center[3*i_face  ];
      cell_center[3*i_cell+1] += face_center[3*i_face+1];
      cell_center[3*i_cell+2] += face_center[3*i_face+2];
    }
    cell_center[3*i_cell  ] = cell_center[3*i_cell  ] * pond;
    cell_center[3*i_cell+1] = cell_center[3*i_cell+1] * pond;
    cell_center[3*i_cell+2] = cell_center[3*i_cell+2] * pond;
  }

  // Dual
  for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

    double xc = cell_center[3*i_cell  ];
    double yc = cell_center[3*i_cell+1];
    double zc = cell_center[3*i_cell+2];

    for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {

      int i_face   = PDM_ABS (cell_face[idx_face])-1;
      int sgn_face = PDM_SIGN(cell_face[idx_face]);
      double xf = face_center[3*i_face  ];
      double yf = face_center[3*i_face+1];
      double zf = face_center[3*i_face+2];

      for(int idx_edge = face_edge_idx[i_face]; idx_edge < face_edge_idx[i_face+1]; ++idx_edge) {

        int i_edge   = PDM_ABS (face_edge[idx_edge])-1;
        int sgn_edge = PDM_SIGN(face_edge[idx_edge]);
        double xe = edge_center[3*i_edge  ];
        double ye = edge_center[3*i_edge+1];
        double ze = edge_center[3*i_edge+2];

        double ux = xf - xc;
        double uy = yf - yc;
        double uz = zf - zc;

        double vx = xe - xc;
        double vy = ye - yc;
        double vz = ze - zc;

        double nx = 0.5 * (uy * vz - uz * vy);
        double ny = 0.5 * (uz * vx - ux * vz);
        double nz = 0.5 * (ux * vy - uy * vx);

        edge_surf[3*i_edge  ] += sgn_face * sgn_edge * nx;
        edge_surf[3*i_edge+1] += sgn_face * sgn_edge * ny;
        edge_surf[3*i_edge+2] += sgn_face * sgn_edge * nz;

        int i_vtx1 = PDM_ABS (edge_vtx[2*i_edge  ])-1;
        int i_vtx2 = PDM_ABS (edge_vtx[2*i_edge+1])-1;

        double dx = coords[3*i_vtx2  ] - coords[3*i_vtx1  ];
        double dy = coords[3*i_vtx2+1] - coords[3*i_vtx1+1];
        double dz = coords[3*i_vtx2+2] - coords[3*i_vtx1+2];

        double test = dx * nx + dy * ny + dz *nz;

        // printf("test = %i %12.5e %i %i \n", i_edge, sgn_face * sgn_edge * test, sgn_face, sgn_edge);
        printf("i_cell = %i i_face = %i i_edge = %i -> %12.5e \n", i_cell, cell_face[idx_face], face_edge[idx_edge], sgn_face * sgn_edge * test);
        // printf("test = %12.5e \n", test);
        // printf("xc, yc, zc = %12.5e %12.5e %12.5e | xf, yf, zf = %12.5e %12.5e %12.5e \n", xc, yc, zc, xf, yf, zf);
        // printf("xe, ye, ze = %12.5e %12.5e %12.5e | xf, yf, zf = %12.5e %12.5e %12.5e \n", xe, ye, ze, xf, yf, zf);
        // printf("ux, uy, uz = %12.5e %12.5e %12.5e | vx, vy, vz = %12.5e %12.5e %12.5e \n", ux, uy, uz, vx, vy, vz);
        // printf("nx, ny, nz = %12.5e %12.5e %12.5e | dx, dy, dz = %12.5e %12.5e %12.5e \n", nx, ny, nz, dx, dy, dz);


      } /* End face_edge */
    } /* End cell_face */
  } /* End cell */

  PDM_log_trace_array_int (face_edge_idx, n_face+1             , "dface_edge_idx ::");
  PDM_log_trace_array_int (face_edge    , face_edge_idx[n_face], "face_edge ::");
  PDM_log_trace_array_int (edge_vtx     , 2* n_edge            , "edge_vtx ::");

  double *vtx_dual_surf = malloc (sizeof(double) * 3 * n_vtx);
  for (int i = 0; i < 3*n_vtx; i++) {
    vtx_dual_surf[i] = 0.;
  }

  for (int i_edge = 0; i_edge < n_edge; i_edge++) {
    int sgn_vtx = 1;
    for (int idx_vtx = 2*i_edge; idx_vtx < 2*(i_edge+1); idx_vtx++) {

      int i_vtx = PDM_ABS(edge_vtx[idx_vtx]) - 1;
      for (int i = 0; i < 3; i++) {
        vtx_dual_surf[3*i_vtx + i] += sgn_vtx * edge_surf[3*i_edge + i];
      }

      sgn_vtx = -sgn_vtx;
    }
  }

  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    double mag = PDM_MODULE (vtx_dual_surf + 3*i_vtx);
    if (mag < 1e-13) {
      printf("vtx %d: (%f %f %f)   mag = %g\n", i_vtx, coords[3*i_vtx], coords[3*i_vtx+1], coords[3*i_vtx+2], mag);
    }
  }
  free(vtx_dual_surf);


  free(edge_surf);
  free(face_center);
  free(cell_center);
  free(edge_center);
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
  PDM_g_num_t nx = 10;
  PDM_g_num_t ny = 10;
  PDM_g_num_t nz = 10;

  double lengthx = 1.;
  double lengthy = 1.;
  double lengthz = 1.;

  int n_part      = 1;
  int post        = 0;
  int randomize   = 0;
  int random_seed = 0;

  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &nx,
              &ny,
              &nz,
              &lengthx,
              &lengthy,
              &lengthz,
              &n_part,
              &randomize,
              &random_seed,
              &post,
      (int *) &method);


  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  /*
   *  Create distributed mesh
   */
  double xmin = 0.;
  double ymin = 0.;
  double zmin = 0.;

  PDM_g_num_t  ng_cell;
  PDM_g_num_t  ng_face;
  PDM_g_num_t  ng_vtx;
  int          n_face_group;
  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int         *dcell_face_idx  = NULL;
  PDM_g_num_t *dcell_face      = NULL;
  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  PDM_poly_vol_gen (comm,
                    xmin,
                    ymin,
                    zmin,
                    lengthx,
                    lengthy,
                    lengthz,
                    nx,
                    ny,
                    nz,
                    randomize,
                    random_seed,
                    &ng_cell,
                    &ng_face,
                    &ng_vtx,
                    &n_face_group,
                    &dn_cell,
                    &dn_face,
                    &dn_vtx,
                    &dcell_face_idx,
                    &dcell_face,
                    &dface_cell,
                    &dface_vtx_idx,
                    &dface_vtx,
                    &dvtx_coord,
                    &dface_group_idx,
                    &dface_group);


  /*
   *  Create mesh partitions
   */
  PDM_g_num_t* cell_distribution = PDM_compute_entity_distribution(comm, dn_cell);
  PDM_g_num_t* face_distribution = PDM_compute_entity_distribution(comm, dn_face);
  PDM_g_num_t* part_distribution = PDM_compute_entity_distribution(comm, n_part );
  PDM_g_num_t* vtx_distribution  = PDM_compute_entity_distribution(comm, dn_vtx );

  /*
   * Generate edge numbering
   */
  int n_edge_elt_tot = dface_vtx_idx[dn_face];
  PDM_g_num_t* tmp_dface_edge         = (PDM_g_num_t *) malloc(     n_edge_elt_tot    * sizeof(PDM_g_num_t) );
  int*         tmp_parent_elmt_pos    = (int         *) malloc(     n_edge_elt_tot    * sizeof(int        ) );
  int*         tmp_dface_edge_vtx_idx = (int         *) malloc( ( n_edge_elt_tot + 1) * sizeof(int        ) );
  PDM_g_num_t* tmp_dface_edge_vtx     = (PDM_g_num_t *) malloc( 2 * n_edge_elt_tot    * sizeof(PDM_g_num_t) );

  int n_elmt_current = 0;
  int n_edge_current = 0;
  tmp_dface_edge_vtx_idx[0] = 0;
  PDM_poly2d_decomposes_edges(dn_face,
                              &n_elmt_current,
                              &n_edge_current,
                              face_distribution[i_rank],
                              -1,
                              dface_vtx,
                              dface_vtx_idx,
                              tmp_dface_edge_vtx_idx,
                              tmp_dface_edge_vtx,
                              tmp_dface_edge,
                              NULL,
                              NULL,
                              tmp_parent_elmt_pos);
  assert(n_edge_current == n_edge_elt_tot);
  free(tmp_parent_elmt_pos);

  int  dn_edge = -1;
  PDM_g_num_t  *dedge_distrib;
  int          *dedge_vtx_idx;
  PDM_g_num_t  *dedge_vtx;
  int          *dedge_face_idx;
  PDM_g_num_t  *dedge_face;

  PDM_generate_entitiy_connectivity_raw(comm,
                                        vtx_distribution[n_rank],
                                        n_edge_elt_tot,
                                        tmp_dface_edge,
                                        tmp_dface_edge_vtx_idx,
                                        tmp_dface_edge_vtx,
                                        &dn_edge,
                                        &dedge_distrib,
                                        &dedge_vtx_idx,
                                        &dedge_vtx,
                                        &dedge_face_idx,
                                        &dedge_face);

  int         *dface_edge_idx = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  PDM_dconnectivity_transpose (comm,
                               dedge_distrib,
                               face_distribution,
                               dedge_face_idx,
                               dedge_face,
                               1,
                               &dface_edge_idx,
                               &dface_edge);


  _compute_dual_closed (dn_cell,
                        dn_face,
                        dn_edge,
                        dn_vtx,
                        dcell_face_idx,
                        dcell_face,
                        dface_edge_idx,
                        dface_edge,
                        dedge_vtx_idx,
                        dedge_vtx,
                        dvtx_coord);




  PDM_MPI_Finalize();

  return 0;
}
