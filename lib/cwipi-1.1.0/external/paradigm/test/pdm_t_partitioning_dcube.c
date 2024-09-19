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
#include "pdm_dmesh_partitioning.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dcube_gen.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
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

static
void
_compute_triangle_surf
(
 double x1,
 double y1,
 double z1,
 double x2,
 double y2,
 double z2,
 double x3,
 double y3,
 double z3,
 double *nx,
 double *ny,
 double *nz
)
{
  double ux = x2 - x1;
  double uy = y2 - y1;
  double uz = z2 - z1;

  double vx = x3 - x1;
  double vy = y3 - y1;
  double vz = z3 - z1;

  *nx = 0.5 * (uy * vz - uz * vy);
  *ny = 0.5 * (uz * vx - ux * vz);
  *nz = 0.5 * (ux * vy - uy * vx);
}


static
void
compute_dual_mesh_metrics
(
 int       n_part,
 int      *pn_vtx,
 int      *pn_cell,
 int      *pn_faces,
 int      *pn_edge,
 double  **pvtx_coord,
 int     **pface_cell,
 int     **pface_vtx_idx,
 int     **pface_vtx,
 int     **pedge_face_idx,
 int     **pedge_face,
 int     **pedge_vtx,
 int       n_face_group,
 int     **pface_group_idx,
 int     **pface_group,
 double ***pedge_surf,
 double ***pdual_vol
)
{
  /*
   * Compute cell_center and face center coordinates
   */
  double** center_cell        = (double **) malloc( n_part * sizeof(double *) );
  double** center_face        = (double **) malloc( n_part * sizeof(double *) );
  double** surface_face       = (double **) malloc( n_part * sizeof(double *) );
  double** ponderate_face_vtx = (double **) malloc( n_part * sizeof(double *) );

  for (int i_part = 0; i_part < n_part; i_part++){

    center_cell [i_part] = (double *) malloc( 3 * pn_cell[i_part]  * sizeof(double) );
    center_face [i_part] = (double *) malloc( 3 * pn_faces[i_part] * sizeof(double) );
    surface_face[i_part] = (double *) malloc( 3 * pn_faces[i_part] * sizeof(double) );

    ponderate_face_vtx[i_part] = (double *) malloc( pface_vtx_idx[i_part][pn_faces[i_part]] * sizeof(double) );
    int* count_cell      = (int    *) malloc(     pn_cell[i_part]  * sizeof(int   ) );

    for (int iface = 0 ; iface < pn_faces[i_part]; iface++) {
      center_face[i_part][3*iface  ] = 0.;
      center_face[i_part][3*iface+1] = 0.;
      center_face[i_part][3*iface+2] = 0.;
    }

    for (int icell = 0 ; icell < pn_cell[i_part]; icell++) {
      center_cell[i_part][3*icell  ] = 0.;
      center_cell[i_part][3*icell+1] = 0.;
      center_cell[i_part][3*icell+2] = 0.;
      count_cell[icell]              = 0;
    }

    for (int iface = 0 ; iface < pn_faces[i_part]; iface++) {
      double x2 = 0.;
      double y2 = 0.;
      double z2 = 0.;
      int n_vtx_on_face = pface_vtx_idx[i_part][iface+1] - pface_vtx_idx[i_part][iface];
      for(int idx_vtx = pface_vtx_idx[i_part][iface]; idx_vtx < pface_vtx_idx[i_part][iface+1]; ++idx_vtx) {
        int i_vtx = PDM_ABS(pface_vtx[i_part][idx_vtx])-1;
        x2 += pvtx_coord[i_part][3*i_vtx  ];
        y2 += pvtx_coord[i_part][3*i_vtx+1];
        z2 += pvtx_coord[i_part][3*i_vtx+2];
      }

      center_face[i_part][3*iface  ] = x2/n_vtx_on_face;
      center_face[i_part][3*iface+1] = y2/n_vtx_on_face;
      center_face[i_part][3*iface+2] = z2/n_vtx_on_face;

      int i_cell1 = PDM_ABS(pface_cell[i_part][2*iface  ])-1;
      int i_cell2 = PDM_ABS(pface_cell[i_part][2*iface+1])-1;

      // printf(" [%i] --> icell1 = %i | icell2 = %i [%i] \n", i_part, i_cell1, i_cell2, pn_cell[i_part]);
      center_cell[i_part][3*i_cell1  ] += x2;
      center_cell[i_part][3*i_cell1+1] += y2;
      center_cell[i_part][3*i_cell1+2] += z2;
      count_cell[i_cell1] += n_vtx_on_face;

      if(i_cell2 > 0 ){
        center_cell[i_part][3*i_cell2  ] += x2;
        center_cell[i_part][3*i_cell2+1] += y2;
        center_cell[i_part][3*i_cell2+2] += z2;
        count_cell[i_cell2] += n_vtx_on_face;
      }

      /* Surface for faces */
      surface_face[i_part][3*iface  ] = 0.;
      surface_face[i_part][3*iface+1] = 0.;
      surface_face[i_part][3*iface+2] = 0.;

      double xcf = center_face[i_part][3*iface  ];
      double ycf = center_face[i_part][3*iface+1];
      double zcf = center_face[i_part][3*iface+2];

      int beg = pface_vtx_idx[i_part][iface];

      /* Init */
      for(int idx_vtx2 = 0; idx_vtx2 < n_vtx_on_face; ++idx_vtx2) {
        ponderate_face_vtx[i_part][beg+idx_vtx2] = 0.;
      }

      for(int idx_vtx2 = 0; idx_vtx2 < n_vtx_on_face; ++idx_vtx2) {
        int idx_vtx3 = (idx_vtx2 + 1) % n_vtx_on_face;
        // printf("[%i] idx_vtx : %i %i \n", iface, idx_vtx2, idx_vtx3);
        int i_vtx2 = PDM_ABS(pface_vtx[i_part][beg+idx_vtx2])-1;
        int i_vtx3 = PDM_ABS(pface_vtx[i_part][beg+idx_vtx3])-1;

        double nx, ny, nz;

        double xvtx2 = pvtx_coord[i_part][3*i_vtx2  ];
        double yvtx2 = pvtx_coord[i_part][3*i_vtx2+1];
        double zvtx2 = pvtx_coord[i_part][3*i_vtx2+2];

        double xvtx3 = pvtx_coord[i_part][3*i_vtx3  ];
        double yvtx3 = pvtx_coord[i_part][3*i_vtx3+1];
        double zvtx3 = pvtx_coord[i_part][3*i_vtx3+2];

        _compute_triangle_surf(xcf, ycf, zcf, xvtx2, yvtx2, zvtx2, xvtx3, yvtx3, zvtx3, &nx, &ny, &nz);

        surface_face[i_part][3*iface  ] += nx;
        surface_face[i_part][3*iface+1] += ny;
        surface_face[i_part][3*iface+2] += nz;

        // middle of current edge
        double mid_x = 0.5 * (xvtx2 + xvtx3);
        double mid_y = 0.5 * (yvtx2 + yvtx3);
        double mid_z = 0.5 * (zvtx2 + zvtx3);

        _compute_triangle_surf(xcf, ycf, zcf, xvtx2, yvtx2, zvtx2, mid_x, mid_y, mid_z, &nx, &ny, &nz);
        ponderate_face_vtx[i_part][beg+idx_vtx2] += sqrt(nx*nx + ny*ny + nz*nz);

        _compute_triangle_surf(xcf, ycf, zcf, mid_x, mid_y, mid_z, xvtx3, yvtx3, zvtx3, &nx, &ny, &nz);
        ponderate_face_vtx[i_part][beg+idx_vtx3] += sqrt(nx*nx + ny*ny + nz*nz);

      }

      // Compute norm
      double sn = sqrt(surface_face[i_part][3*iface  ]*surface_face[i_part][3*iface  ] +
                       surface_face[i_part][3*iface+1]*surface_face[i_part][3*iface+1] +
                       surface_face[i_part][3*iface+2]*surface_face[i_part][3*iface+2]);

      double invsn = 1./sn;

      for(int idx_vtx2 = 0; idx_vtx2 < n_vtx_on_face; ++idx_vtx2) {
        ponderate_face_vtx[i_part][beg+idx_vtx2] = ponderate_face_vtx[i_part][beg+idx_vtx2] * invsn;
        printf("ponderate_face_vtx[%i] = %12.5e \n", beg+idx_vtx2, ponderate_face_vtx[i_part][beg+idx_vtx2]);
      }

    }

    /* Adim cell center */
    for (int icell = 0 ; icell < pn_cell[i_part]; icell++) {
      double inv = 1./count_cell[icell];
      center_cell[i_part][3*icell  ] = center_cell[i_part][3*icell  ] * inv;
      center_cell[i_part][3*icell+1] = center_cell[i_part][3*icell+1] * inv;
      center_cell[i_part][3*icell+2] = center_cell[i_part][3*icell+2] * inv;
    }

    for (int iface = 0 ; iface < pn_faces[i_part]; iface++) {
      printf(" center_face[%i] = %12.5e %12.5e %12.5e \n", iface, center_face[i_part][3*iface  ], center_face[i_part][3*iface+1], center_face[i_part][3*iface+2]);
    }
    for (int iface = 0 ; iface < pn_faces[i_part]; iface++) {
      printf(" surface_face[%i] = %12.5e %12.5e %12.5e \n", iface, surface_face[i_part][3*iface  ], surface_face[i_part][3*iface+1], surface_face[i_part][3*iface+2]);
    }
    for (int icell = 0 ; icell < pn_cell[i_part]; icell++) {
      printf(" center_cell[%i] = %12.5e %12.5e %12.5e (%i) \n", icell, center_cell[i_part][3*icell  ], center_cell[i_part][3*icell+1], center_cell[i_part][3*icell+2], count_cell[icell]);
    }

    free(count_cell);
  }

  /*
   * Compute normal associate to edge
   */
  double** edge_surf = (double **) malloc( n_part * sizeof(double *) );
  double** dual_vol  = (double **) malloc( n_part * sizeof(double *) );
  for (int i_part = 0; i_part < n_part; i_part++){

    int    *_pedge_vtx  = pedge_vtx [i_part];
    double *_pvtx_coord = pvtx_coord[i_part];

    edge_surf[i_part]  = (double *) malloc( 3 * pn_edge[i_part] * sizeof(double *) );
    dual_vol [i_part]  = (double *) malloc(     pn_vtx [i_part] * sizeof(double *) );
    double *_edge_surf = edge_surf[i_part];
    double *_dual_vol  = dual_vol[i_part];

    for(int ivtx = 0; ivtx < pn_vtx[i_part]; ++ivtx) {
      _dual_vol[ivtx] = 0.;
    }

    for (int iedge=0 ; iedge < pn_edge[i_part]; iedge++) {

      // Chaque edge est lié a deux vertex -_-
      int i_vtx1 = _pedge_vtx[2*iedge  ]-1;
      int i_vtx2 = _pedge_vtx[2*iedge+1]-1;

      _edge_surf[3*iedge  ] = 0.;
      _edge_surf[3*iedge+1] = 0.;
      _edge_surf[3*iedge+2] = 0.;

      double x1 = 0.5 * ( _pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
      double y1 = 0.5 * ( _pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
      double z1 = 0.5 * ( _pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);

      // printf(" -------------------------------------------------------- \n");
      // printf(" x1 = %12.5e | y2 = %12.5e | z2 = %12.5e \n", x1, y1, z1);

      // printf(" edge %i -> %i %i -> n_face = %i \n", iedge, i_vtx1+1, i_vtx2+1, pedge_face_idx[i_part][iedge+1]-pedge_face_idx[i_part][iedge]);

      for(int idx_face = pedge_face_idx[i_part][iedge]; idx_face < pedge_face_idx[i_part][iedge+1]; ++idx_face ) {

        int i_face = PDM_ABS (pedge_face[i_part][idx_face])-1;
        int sgn    = PDM_SIGN(pedge_face[i_part][idx_face]);
        double nx, ny, nz;

        int i_cell1 = PDM_ABS(pface_cell[i_part][2*i_face  ])-1;
        int i_cell2 = PDM_ABS(pface_cell[i_part][2*i_face+1])-1;

        // printf(" i_cell1 = %i | i_cell2 = %i | %i / %i\n", i_cell1, i_cell2, pface_cell[i_part][2*i_face  ], pface_cell[i_part][2*i_face+1]);
        // for(int itmp = pcell_face_idx[i_part][i_cell1]; itmp < pcell_face_idx[i_part][i_cell1+1]; ++itmp) {
        //   printf(" cell_face[%i] = %i \n", itmp, pcell_face[i_part][itmp]);
        // }
        // printf("face_vtx[%i] = ", i_face);
        // for(int itmp = pface_vtx_idx[i_part][i_face]; itmp < pface_vtx_idx[i_part][i_face+1]; ++itmp) {
        //   printf("%i ", pface_vtx[i_part][itmp]);
        // }
        // printf("\n");

        double x2, y2, z2;
        double x3, y3, z3;

        x2 = center_face[i_part][3*i_face  ];
        y2 = center_face[i_part][3*i_face+1];
        z2 = center_face[i_part][3*i_face+2];

        x3 = center_cell[i_part][3*i_cell1  ];
        y3 = center_cell[i_part][3*i_cell1+1];
        z3 = center_cell[i_part][3*i_cell1+2];

        _compute_triangle_surf(x1, y1, z1, x2, y2, z2, x3, y3, z3, &nx, &ny, &nz);
        _edge_surf[3*iedge  ] -= sgn * nx;
        _edge_surf[3*iedge+1] -= sgn * ny;
        _edge_surf[3*iedge+2] -= sgn * nz;

        // printf(" icell1 : sgn = %i | nx = %12.5e | ny = %12.5e | nz = %12.5e\n", sgn, nx, ny, nz);

        double tri_face_center[3];
        tri_face_center[0] = 1./3. * ( x1 + x2 + x3);
        tri_face_center[1] = 1./3. * ( y1 + y2 + y3);
        tri_face_center[2] = 1./3. * ( z1 + z2 + z3);

        double dh_vtx1[3];
        dh_vtx1[0] = _pvtx_coord[3*i_vtx1  ] - tri_face_center[0];
        dh_vtx1[1] = _pvtx_coord[3*i_vtx1+1] - tri_face_center[1];
        dh_vtx1[2] = _pvtx_coord[3*i_vtx1+2] - tri_face_center[2];

        double dh_vtx2[3];
        dh_vtx2[0] = _pvtx_coord[3*i_vtx2  ] - tri_face_center[0];
        dh_vtx2[1] = _pvtx_coord[3*i_vtx2+1] - tri_face_center[1];
        dh_vtx2[2] = _pvtx_coord[3*i_vtx2+2] - tri_face_center[2];

        double surface_vector[3];
        surface_vector[0] = nx;
        surface_vector[1] = ny;
        surface_vector[2] = nz;

        double dual_vol1 = 1./3. * PDM_DOT_PRODUCT(dh_vtx1, surface_vector);
        double dual_vol2 = 1./3. * PDM_DOT_PRODUCT(dh_vtx2, surface_vector);

        _dual_vol[i_vtx1] += PDM_ABS(dual_vol1);
        _dual_vol[i_vtx2] += PDM_ABS(dual_vol2);

        if(i_cell2 > -1 ){

          x3 = center_cell[i_part][3*i_cell2  ];
          y3 = center_cell[i_part][3*i_cell2+1];
          z3 = center_cell[i_part][3*i_cell2+2];


          _compute_triangle_surf(x1, y1, z1, x2, y2, z2, x3, y3, z3, &nx, &ny, &nz);

          // printf(" icell2 : sgn = %i | nx = %12.5e | ny = %12.5e | nz = %12.5e\n", sgn, nx, ny, nz);

          _edge_surf[3*iedge  ] += sgn * nx;
          _edge_surf[3*iedge+1] += sgn * ny;
          _edge_surf[3*iedge+2] += sgn * nz;

          tri_face_center[0] = 1./3. * ( x1 + x2 + x3);
          tri_face_center[1] = 1./3. * ( y1 + y2 + y3);
          tri_face_center[2] = 1./3. * ( z1 + z2 + z3);

          surface_vector[0] = nx;
          surface_vector[1] = ny;
          surface_vector[2] = nz;

          double dual_vol1_cell2 = 1./3. * PDM_DOT_PRODUCT(dh_vtx1, surface_vector);
          double dual_vol2_cell2 = 1./3. * PDM_DOT_PRODUCT(dh_vtx2, surface_vector);

          _dual_vol[i_vtx1] += PDM_ABS(dual_vol1_cell2);
          _dual_vol[i_vtx2] += PDM_ABS(dual_vol2_cell2);

        }

        // double surf_norm = sqrt(_edge_surf[3*iedge  ]*_edge_surf[3*iedge  ] +  _edge_surf[3*iedge+1]*_edge_surf[3*iedge+1] + _edge_surf[3*iedge+2]*_edge_surf[3*iedge+2]);
        // printf(" edge_surf[%i] = %12.5e %12.5e %12.5e -->  %12.5e (%i) \n", iedge, _edge_surf[3*iedge  ], _edge_surf[3*iedge+1], _edge_surf[3*iedge+2], surf_norm, pedge_face_idx[i_part][iedge+1]-pedge_face_idx[i_part][iedge]);

      }

      // double surf_norm = sqrt(_edge_surf[3*iedge  ]*_edge_surf[3*iedge  ] +  _edge_surf[3*iedge+1]*_edge_surf[3*iedge+1] + _edge_surf[3*iedge+2]*_edge_surf[3*iedge+2]);
      // printf(" edge_surf[%i] = %12.5e %12.5e %12.5e -->  %12.5e (%i) \n", iedge, _edge_surf[3*iedge  ], _edge_surf[3*iedge+1], _edge_surf[3*iedge+2], surf_norm, pedge_face_idx[i_part][iedge+1]-pedge_face_idx[i_part][iedge]);

    }

    for (int iedge=0 ; iedge < pn_edge[i_part]; iedge++) {
      double surf_norm = sqrt(_edge_surf[3*iedge  ]*_edge_surf[3*iedge  ] +  _edge_surf[3*iedge+1]*_edge_surf[3*iedge+1] + _edge_surf[3*iedge+2]*_edge_surf[3*iedge+2]);
      assert(surf_norm > 0.);
      printf(" edge_surf[%i] = %12.5e %12.5e %12.5e -->  %12.5e (%i) \n", iedge, _edge_surf[3*iedge  ], _edge_surf[3*iedge+1], _edge_surf[3*iedge+2], surf_norm, pedge_face_idx[i_part][iedge+1]-pedge_face_idx[i_part][iedge]);
    }

    double tot_volume = 0.;
    for(int ivtx = 0; ivtx < pn_vtx[i_part]; ++ivtx) {
      printf(" _dual_vol[%i] = %12.5e \n", ivtx, _dual_vol[ivtx]);
      tot_volume += _dual_vol[ivtx];
    }
    printf(" tot_volume = %12.5e \n", tot_volume);

    /* Flux balance */
    double *flux_bal = (double *) malloc( 3 * pn_vtx [i_part] * sizeof(double *) );
    for(int ivtx = 0; ivtx < pn_vtx[i_part]; ++ivtx) {
      flux_bal[3*ivtx  ] = 0.;
      flux_bal[3*ivtx+1] = 0.;
      flux_bal[3*ivtx+2] = 0.;
    }

    /* Manage all interior */
    for (int iedge=0 ; iedge < pn_edge[i_part]; iedge++) {

      int i_vtx1 = _pedge_vtx[2*iedge  ]-1;
      int i_vtx2 = _pedge_vtx[2*iedge+1]-1;

      printf("_edge_surf[%i] = %12.5e -> flux_bal[%i] = %12.5e |  flux_bal[%i] = %12.5e | \n",
             iedge, _edge_surf[3*iedge  ], i_vtx1, flux_bal[3*i_vtx1  ], i_vtx2, flux_bal[3*i_vtx2  ]);

      flux_bal[3*i_vtx1  ] += _edge_surf[3*iedge  ];
      flux_bal[3*i_vtx1+1] += _edge_surf[3*iedge+1];
      flux_bal[3*i_vtx1+2] += _edge_surf[3*iedge+2];

      flux_bal[3*i_vtx2  ] -= _edge_surf[3*iedge  ];
      flux_bal[3*i_vtx2+1] -= _edge_surf[3*iedge+1];
      flux_bal[3*i_vtx2+2] -= _edge_surf[3*iedge+2];
    }

    /* Manage Bnd */
    for(int i_group = 0; i_group < n_face_group; ++i_group) {
      for(int idx_face = pface_group_idx[i_part][i_group]; idx_face < pface_group_idx[i_part][i_group+1]; ++idx_face ) {

        int iface = pface_group[i_part][idx_face]-1;

        int n_vtx_on_face = pface_vtx_idx[i_part][iface+1] - pface_vtx_idx[i_part][iface];
        double pond = 1./n_vtx_on_face;
        for(int idx_vtx = pface_vtx_idx[i_part][iface]; idx_vtx < pface_vtx_idx[i_part][iface+1]; ++idx_vtx) {
          int i_vtx = pface_vtx[i_part][idx_vtx]-1;
          flux_bal[3*i_vtx  ] += pond * surface_face[i_part][3*iface  ];
          flux_bal[3*i_vtx+1] += pond * surface_face[i_part][3*iface+1];
          flux_bal[3*i_vtx+2] += pond * surface_face[i_part][3*iface+2];
        }

      }
    }


    for(int ivtx = 0; ivtx < pn_vtx[i_part]; ++ivtx) {
      double surf_bal_norm = sqrt(flux_bal[3*ivtx  ]*flux_bal[3*ivtx  ] +  flux_bal[3*ivtx+1]*flux_bal[3*ivtx+1] + flux_bal[3*ivtx+2]*flux_bal[3*ivtx+2]);
      printf(" flux_bal[%i] = %12.5e %12.5e %12.5e -->  %12.5e \n", ivtx, flux_bal[3*ivtx  ], flux_bal[3*ivtx+1], flux_bal[3*ivtx+2], surf_bal_norm);
    }

    free(flux_bal);

  }

  for (int i_part=0; i_part < n_part; i_part++){
    free(center_cell[i_part]);
    free(center_face[i_part]);
    free(surface_face[i_part]);
    free(ponderate_face_vtx[i_part]);
  }
  free(center_cell);
  free(center_face);
  free(surface_face);
  free(ponderate_face_vtx);

  *pedge_surf = edge_surf;
  *pdual_vol  = dual_vol;
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
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PARMETIS;

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

  struct timeval t_elaps_debut;

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

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

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;

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

  if (0 == 1) {

    PDM_printf("[%i] n_face_group    : %i\n", i_rank, n_face_group);
    PDM_printf("[%i] dn_cell        : %i\n", i_rank, dn_cell);
    PDM_printf("[%i] dn_face        : %i\n", i_rank, dn_face);
    PDM_printf("[%i] dn_vtx         : %i\n", i_rank, dn_vtx);

    PDM_printf("[%i] dface_cell     : ", i_rank);
    for (int i = 0; i < 2 * dn_face; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_cell[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_vtx_idx   : ", i_rank);
    for (int i = 0; i < dn_face + 1; i++)
      PDM_printf(" %i", dface_vtx_idx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_vtx      : ", i_rank);
    for (int i = 0; i < dface_vtx_idx[dn_face]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_vtx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dvtx_coord     : ", i_rank);
    for (int i = 0; i < 3*dn_vtx; i++)
      PDM_printf(" %12.5e", dvtx_coord[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_group_idx : ", i_rank);
    for (int i = 0; i < n_face_group + 1; i++)
      PDM_printf(" %i", dface_group_idx[i]);
    PDM_printf("\n");

    PDM_printf("[%i] dface_group    : ", i_rank);
    for (int i = 0; i < dface_group_idx[n_face_group]; i++)
      PDM_printf(" "PDM_FMT_G_NUM, dface_group[i]);
    PDM_printf("\n");

  }

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

  /* Make ascending connectivity */
  int          *dface_edge_idx;
  PDM_g_num_t  *dface_edge;
  if (post) {
    PDM_log_trace_array_long(face_distribution, n_rank+1, "face_distribution::");
    PDM_log_trace_array_long(dedge_distrib, n_rank+1, "dedge_distrib::");
    PDM_log_trace_array_int(dedge_vtx_idx, dn_edge+1, "dedge_vtx_idx::");
    PDM_log_trace_array_long(dedge_vtx, dedge_vtx_idx[dn_edge], "dedge_vtx::");
    PDM_log_trace_array_int(dedge_face_idx, dn_edge, "dedge_face_idx::");
    PDM_log_trace_array_long(dedge_face, dedge_face_idx[dn_edge], "dedge_face::");
  }

  PDM_dconnectivity_transpose(comm,
                              dedge_distrib,
                              face_distribution,
                              dedge_face_idx,
                              dedge_face,
                              1,
                              &dface_edge_idx,
                              &dface_edge);

  if (post) {
    printf("dn_edge = %i \n", dn_edge);
    PDM_log_trace_array_int(dface_edge_idx, dn_face, "dface_edge_idx::");
    PDM_log_trace_array_long(dface_edge, dface_edge_idx[dn_face], "dface_edge::");
  }

  // int flags = PDM_PART_FACE_CELL|PDM_PART_CELL_FACE;
  // printf("PDM_HASFLAG(flags, PDM_PART_FACE_CELL) :: %d\n", PDM_HASFLAG(flags, PDM_PART_FACE_CELL) );
  // printf("PDM_HASFLAG(flags, PDM_PART_CELL_FACE) :: %d\n", PDM_HASFLAG(flags, PDM_PART_CELL_FACE) );
  // printf("PDM_HASFLAG(flags, PDM_PART_FACE_VTX) :: %d\n" , PDM_HASFLAG(flags, PDM_PART_FACE_VTX) );
  // printf("x::PDM_HASFLAG(flags, PDM_PART_FACE_VTX) :: %x\n", PDM_PART_FACE_VTX);

  gettimeofday(&t_elaps_debut, NULL);

  // printf("part_distribution::\n");
  // for(int i_part = 0; i_part < n_rank+1; ++i_part){
  //   printf("%d ", part_distribution[i_part]);
  // }
  // printf("\n");

  /*
   * Compute dual graph
   */
  PDM_g_num_t* dual_graph_idx;
  PDM_g_num_t* dual_graph;
  int* dcell_face_idx;
  PDM_g_num_t* dcell_face;
  PDM_para_graph_dual_from_arc2node(comm,
                                    cell_distribution,
                                    face_distribution,
                                    dface_cell,
                    (PDM_g_num_t**) &dual_graph_idx,
                    (PDM_g_num_t**) &dual_graph,
                                    PDM_SPLIT_DUAL_WITH_PARMETIS,
                    (int        **) &dcell_face_idx,
                    (PDM_g_num_t**) &dcell_face);

  // Test graph creation from cell_face connectivity
  if (0 == 1) {

    if (0 == 1) {
      printf("dcell_face_idx :");
      for (int i =0; i < dn_cell+1; i++)
        printf(" %d", dcell_face_idx[i]);
      printf("\n");

      printf("dcell_face:: %d \n", dn_cell);
      for(int i = 0; i < dn_cell; ++i){
        printf("Local cell %d :", i);
        for(int j = dcell_face_idx[i]; j < dcell_face_idx[i+1]; ++j){
          printf(" "PDM_FMT_G_NUM"", dcell_face[j]);
        }
        PDM_printf("\n");
      }
    }

    free(dual_graph_idx);
    free(dual_graph);
    dual_graph_idx = NULL;
    dual_graph = NULL;

    PDM_para_graph_dual_from_node2arc(comm,
                                      cell_distribution,
                                      face_distribution,
                                      dcell_face_idx,
                                      dcell_face,
                      (PDM_g_num_t**) &dual_graph_idx,
                      (PDM_g_num_t**) &dual_graph);
  }

  if (post) {
    PDM_log_trace_array_long(dual_graph_idx, dn_cell+1              , "pdm_t_partitioning_dcube::dual_graph_idx::");
    PDM_log_trace_array_long(dual_graph    , dual_graph_idx[dn_cell], "pdm_t_partitioning_dcube::dual_graph::");
  }

  // free(dual_graph_idx);
  // free(dual_graph);
  // mpirun -np 2 ./paradigm/test/pdm_t_partitioning_dcube -n 23 -n_part 1 -parmetis
  // PDM_para_graph_dual_from_combine_connectivity(comm,
  //                                               cell_distribution,
  //                                               face_distribution,
  //                                               vtx_distribution,
  //                                               dcell_face_idx,
  //                                               dcell_face,
  //                                               dface_vtx_idx,
  //                                               dface_vtx,
  //                              (PDM_g_num_t**) &dual_graph_idx,
  //                              (PDM_g_num_t**) &dual_graph);

  /*
   * Split it !!! CAUTION dn_cell can be different of the size of dual graph !!!
   */
  // printf("PDM_split_graph\n");
  int* cell_part    = (int *) malloc( sizeof(int) * dn_cell );
  int* dcell_weight = (int *) malloc( sizeof(int) * dn_cell );
  for(int i = 0; i < dn_cell; ++i){
    dcell_weight[i] = dual_graph_idx[i+1] - dual_graph_idx[i];
  }

  if( 0 == 1 ){
    printf("n_cell_block:: %d \n", dn_cell);
    for(int i = 0; i < dn_cell; ++i){
      printf(" dual_graph_idx = "PDM_FMT_G_NUM" ---> \n", dual_graph_idx[i]);
      for(int i_data = dual_graph_idx[i]; i_data < dual_graph_idx[i+1]; ++i_data){
        // printf("%d ", dual_graph[i_data]);
        printf("\t dual_graph[%d] = "PDM_FMT_G_NUM" \n", i_data, dual_graph[i_data]);
      }
      printf("\n");
    }
  }

  int tn_part = part_distribution[n_rank];

  double *part_frac = NULL;
  if (1 == 0) {
    part_frac = (double *) malloc(sizeof(double) * tn_part );
    for (int i_part = 0; i_part < tn_part-1; i_part++)
    {
      if (i_part % 2 == 0) part_frac[i_part] = (double) 0.5*(1./tn_part);
      else                 part_frac[i_part] = (double) 1.5*(1./tn_part);
    }
    if (tn_part % 2 == 0) part_frac[tn_part-1] = (double) 1.5*(1./tn_part);
    else                  part_frac[tn_part-1] = (double) 1.*(1./tn_part);
    PDM_printf("Testing with heterogeneous part sizes\n");
  }

  PDM_para_graph_split (part_method,
                        cell_distribution,
                        dual_graph_idx,
                        dual_graph,
                        NULL,
                        NULL,
                        tn_part,
                        part_frac, // Or NULL for homogeneous parts
                        cell_part,
                        comm);

  // abort();
  if (0 == 1){
    printf("cell_part[%d]::", dn_cell);
    for(int i = 0; i < dn_cell; ++i){
      printf("%d ", cell_part[i]);
    }
    printf("\n");
  }
  if (part_frac != NULL){
    free(part_frac);
  }

  /*
   * On dispose pour chaque cellule de la partition associé : il faut retrouver le
   *  cell_ln_to_gn
   *    Idée : faire un part_to_block avec ln_to_gn == cell_part et distribution imposé = dpart_proc
   *    On devrait recupérer n_part block de données avec le ln_to_gn !
   */
  PDM_g_num_t** pcell_ln_to_gn;
  int*  pn_cell;

  int n_res_part = PDM_part_assemble_partitions(comm,
                                                part_distribution,
                                                cell_distribution,
                                                cell_part,
                                                NULL,
                                                NULL,
                                    (int ** )  &pn_cell,
                            (PDM_g_num_t ***)  &pcell_ln_to_gn,
                                                NULL);
  /*
   * Tentative extented partition
   */
  // PDM_g_num_t** pcell_ln_to_gn_extented;
  // int*          pn_cell_extented;

  // PDM_extend_mesh(comm,
  //                 part_distribution,
  //                 cell_distribution,
  //                 cell_part,
  //                 n_res_part,
  //                 dual_graph_idx,
  //                 dual_graph,
  //                 pn_cell,
  //                 pcell_ln_to_gn,
  //                &pn_cell_extented,
  //                &pcell_ln_to_gn_extented);

  // pn_cell        = pn_cell_extented;
  // pcell_ln_to_gn = pcell_ln_to_gn_extented;
  /*
   *  At this stage we have the cell_ln_to_gn :
   *      --> We need to deduce the other if needed
   *           --> face_ln_to_gn
   *           --> vtx_ln_to_gn
   *  We can do by the relationship of each varibles (connectivity)
   *  For example face_ln_to_gn can be deduce with cell_face connectivity if we have cell_ln_to_gn
   */
  if (post) {
    printf(" pn_cell[0] : %i \n", pn_cell[0]);
    PDM_log_trace_array_long(pcell_ln_to_gn[0]    ,  pn_cell[0]    , "pcell_ln_to_gn::");
  }
  int** pcell_face_idx;
  int** pcell_face;
  int*  pn_faces;
  PDM_g_num_t** pface_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               cell_distribution,
                                               dcell_face_idx,
                                               dcell_face,
                                               n_res_part,
                                               pn_cell,
                      (const PDM_g_num_t ** )  pcell_ln_to_gn,
                            (int         ** ) &pn_faces,
                            (PDM_g_num_t ***) &pface_ln_to_gn,
                            (int         ***) &pcell_face_idx,
                            (int         ***) &pcell_face);

  int **pface_cell;
  PDM_part_reverse_pcellface(n_res_part,
                             pn_cell,
                             pn_faces,
            (const int ** )  pcell_face_idx,
            (const int ** )  pcell_face,
            (      int ***) &pface_cell);

  // int* face_cell_idx = (int *) malloc( (pn_faces[0] + 1 ) * sizeof(int));
  // int* face_cell     = (int *) malloc( (2 * pn_faces[0] ) * sizeof(int));
  // int idx = 0;
  // face_cell_idx[0] = 0;
  // for(int i_face = 0; i_face < pn_faces[0]; ++i_face) {
  //   face_cell_idx[i_face+1] = face_cell_idx[i_face];
  //   if(pface_cell[0][2*i_face + 1 ] == 0) {
  //     face_cell_idx[i_face+1]++;
  //     face_cell[idx++] = pface_cell[0][2*i_face];
  //   } else {
  //     face_cell_idx[i_face+1] += 2;
  //     face_cell[idx++] = pface_cell[0][2*i_face  ];
  //     face_cell[idx++] = pface_cell[0][2*i_face+1];
  //   }
  // }

  // PDM_log_trace_array_int(pface_cell[0]    ,  2 * pn_faces[0]    , "pface_cell::");
  // PDM_log_trace_array_int(pcell_face_idx[0], pn_cell[0]+1                 , "pcell_face_idx::");
  // PDM_log_trace_array_int(pcell_face[0]    , pcell_face_idx[0][pn_cell[0]], "pcell_face::");

  // printf("pn_faces[0] = % i  \n", pn_faces[0]);
  // printf("pn_cell[0]  = % i  \n",pn_cell[0] );
  // PDM_log_trace_array_int(face_cell_idx, pn_faces[0]+1, "face_cell_idx::");
  // PDM_log_trace_array_int(face_cell, face_cell_idx[pn_faces[0]], "face_cell::");
  // free(face_cell_idx);
  // free(face_cell);

  /*
   * Generate vtx
   */

  int** pface_vtx_idx;
  int** pface_vtx;
  int*  pn_vtx;
  PDM_g_num_t** pvtx_ln_to_gn = NULL;

  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               face_distribution,
                                               dface_vtx_idx,
                                               dface_vtx,
                                               n_res_part,
                                               pn_faces,
                        (const PDM_g_num_t **) pface_ln_to_gn,
                           (int         ** )  &pn_vtx,
                           (PDM_g_num_t ***)  &pvtx_ln_to_gn,
                           (int         ***)  &pface_vtx_idx,
                           (int         ***)  &pface_vtx);


  int** pface_edge_idx;
  int** pface_edge;
  int*  pn_edge;
  PDM_g_num_t** pedge_ln_to_gn = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               face_distribution,
                                               dface_edge_idx,
                                               dface_edge,
                                               n_res_part,
                                               pn_faces,
                        (const PDM_g_num_t **) pface_ln_to_gn,
                           (int         ** )  &pn_edge,
                           (PDM_g_num_t ***)  &pedge_ln_to_gn,
                           (int         ***)  &pface_edge_idx,
                           (int         ***)  &pface_edge);

  // Attention il faut modifier le edge face aussi !!!!
  PDM_part_reorient_bound_faces(n_part,
                                pn_faces,
                                pface_cell,
                 (const int **) pcell_face_idx,
                                pcell_face,
                 (const int **) pface_vtx_idx,
                                pface_vtx,
                                pface_edge_idx,
                                pface_edge);

  if (0 == 1){
    for (int i_part=0; i_part < n_res_part; i_part++){
      PDM_printf("[%i] generated facecell part %i [%i]:", i_rank, i_part, pn_faces[i_part]);
      for (int iface=0 ; iface < pn_faces[i_part]; iface++)
        PDM_printf(" %d %d", pface_cell[i_part][2*iface], pface_cell[i_part][2*iface+1]);
      PDM_printf("\n");
    }
  }

  double **pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        vtx_distribution,
                                        dvtx_coord,
                                        pn_vtx,
                (const PDM_g_num_t **)  pvtx_ln_to_gn,
                          (double ***) &pvtx_coord);

  int** pedge_vtx_idx;
  int** pedge_vtx;
  int*  pn_vtx2;
  PDM_g_num_t** pvtx_ln_to_gn2 = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                               dedge_distrib,
                                               dedge_vtx_idx,
                                               dedge_vtx,
                                               n_res_part,
                                               pn_edge,
                      (const PDM_g_num_t ** )  pedge_ln_to_gn,
                            (int         ** ) &pn_vtx2,
                            (PDM_g_num_t ***) &pvtx_ln_to_gn2,
                            (int         ***) &pedge_vtx_idx,
                            (int         ***) &pedge_vtx);
  free(pn_vtx2);
  for(int i_part = 0; i_part < n_res_part; ++i_part) {
    free(pvtx_ln_to_gn2[i_part]);
  }
  free(pvtx_ln_to_gn2);

  int** pedge_face_idx;
  int** pedge_face;
  PDM_part_connectivity_transpose(n_res_part,
                                  pn_faces,
                                  pn_edge,
                                  pface_edge_idx,
                                  pface_edge,
                                  &pedge_face_idx,
                                  &pedge_face);

  if (0 == 1){
    for (int i_part=0; i_part < n_res_part; i_part++){
      PDM_printf("[%i] generated edge_face part %i [%i]: \n", i_rank, i_part, pn_edge[i_part]);
      for (int iedge=0 ; iedge < pn_edge[i_part]; iedge++) {
        PDM_printf(" [%i] -> ", iedge);
        for( int idx_face = pedge_face_idx[i_part][iedge]; idx_face < pedge_face_idx[i_part][iedge+1]; ++idx_face ) {
          PDM_printf(" %i ", pedge_face[i_part][idx_face]);
        }
        PDM_printf("\n");
      }
    }
  }


  /*
   *  Boundary condition (face group )
   */
  PDM_g_num_t** pface_group_ln_to_gn;
  int** pface_group;
  int** pface_group_idx;

  PDM_part_distgroup_to_partgroup(comm,
                                  face_distribution,
                                  n_face_group,
                                  dface_group_idx,
                                  dface_group,
                                  n_res_part,
                                  pn_faces,
           (const PDM_g_num_t **) pface_ln_to_gn,
                (int         ***) &pface_group_idx,
                (int         ***) &pface_group,
                (PDM_g_num_t ***) &pface_group_ln_to_gn);

  PDM_part_reorient_bound_faces(n_part,
                                pn_faces,
                                pface_cell,
                 (const int **) pcell_face_idx,
                      ( int **) pcell_face,
                 (const int **) pface_vtx_idx,
                       (int **) pface_vtx,
                                pface_edge_idx,
                                pface_edge);

  double **edge_surf;
  double **dual_vol;
  compute_dual_mesh_metrics(n_part,
                            pn_vtx,
                            pn_cell,
                            pn_faces,
                            pn_edge,
                            pvtx_coord,
                            pface_cell,
                            pface_vtx_idx,
                            pface_vtx,
                            pedge_face_idx,
                            pedge_face,
                            pedge_vtx,
                            n_face_group,
                            pface_group_idx,
                            pface_group,
                           &edge_surf,
                           &dual_vol);

  for (int i_part=0; i_part < n_res_part; i_part++){
    free(pedge_face_idx[i_part]);
    free(pedge_face[i_part]);
    free(pedge_vtx_idx[i_part]);
    free(pedge_vtx[i_part]);
    free(edge_surf[i_part]);
    free(dual_vol[i_part]);
  }
  free(pedge_face_idx);
  free(pedge_face);
  free(pedge_vtx_idx);
  free(pedge_vtx);
  free(edge_surf);
  free(dual_vol);

  /*
   * Graph communication build
   */
  int** pproc_face_bound_idx;
  int** ppart_face_bound_idx;
  int** pface_bound;

  int **face_is_bnd = (int **) malloc(n_part * sizeof(int*));
  for (int i_part = 0; i_part < n_res_part; i_part++) {
    face_is_bnd[i_part] = (int *) malloc(pn_faces[i_part]*sizeof(int));
    for (int i_face = 0; i_face < pn_faces[i_part]; i_face++){
      if (pface_cell[i_part][2*i_face+1] > 0)
        face_is_bnd[i_part][i_face] = 0;
      else
        face_is_bnd[i_part][i_face] = 1;
    }
  }

  PDM_part_generate_entity_graph_comm(comm,
                                      part_distribution,
                                      face_distribution,
                                      n_part,
                                      pn_faces,
              (const PDM_g_num_t **)  pface_ln_to_gn,
              (const int **)          face_is_bnd,
                           (int ***) &pproc_face_bound_idx,
                           (int ***) &ppart_face_bound_idx,
                           (int ***) &pface_bound,
                                      NULL);
  for (int i_part = 0; i_part < n_res_part; i_part++){
    free(face_is_bnd[i_part]);
  }
  free(face_is_bnd);


  // Attention on veut garder l'orientation donc il y a un signe dans le face_cell / cell_face
  // Reflechir sur les connectivité d'edge également ...

  /*
   * Free
   */
  free(dual_graph_idx);
  free(dual_graph);
  free(cell_part);
  free(dcell_face);
  free(dcell_face_idx);
  free(dcell_weight);
  free(cell_distribution);
  free(face_distribution);
  free(part_distribution);
  free(vtx_distribution);
  for(int i_part = 0; i_part < n_res_part; ++i_part){
    free(pface_ln_to_gn[i_part]);
    free(pcell_ln_to_gn[i_part]);
    // free(pcell_ln_to_gn_extented[i_part]);
    free(pvtx_ln_to_gn[i_part]);
    free(pcell_face[i_part]);
    free(pcell_face_idx[i_part]);
    free(pface_vtx_idx[i_part]);
    free(pface_vtx[i_part]);
    free(pface_group_ln_to_gn[i_part]);
    free(pface_group[i_part]);
    free(pface_group_idx[i_part]);
    free(pproc_face_bound_idx[i_part]);
    free(ppart_face_bound_idx[i_part]);
    free(pface_bound[i_part]);
    free(pface_cell[i_part]);
    free(pvtx_coord[i_part]);
    free(pface_edge_idx[i_part]);
    free(pface_edge[i_part]);
    free(pedge_ln_to_gn[i_part]);
  }
  free(pface_edge_idx);
  free(pface_edge);
  free(pn_edge);
  free(pedge_ln_to_gn);
  free(pcell_face);
  free(pcell_face_idx);
  free(pvtx_ln_to_gn);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pcell_ln_to_gn);
  // free(pcell_ln_to_gn_extented);
  free(pproc_face_bound_idx);
  free(ppart_face_bound_idx);
  free(pface_bound);
  free(pface_ln_to_gn);
  free(pn_cell);
  // free(pn_cell_extented);
  free(pn_faces);
  free(pn_vtx);
  free(pface_group_ln_to_gn);
  free(pface_group);
  free(pface_group_idx);
  free(pface_cell);
  free(pvtx_coord);

  free(dedge_distrib);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face_idx);
  free(dedge_face);

  free(dface_edge_idx);
  free(dface_edge);


  PDM_dcube_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
