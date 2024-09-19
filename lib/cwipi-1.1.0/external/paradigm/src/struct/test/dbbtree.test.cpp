#include "doctest/extensions/doctest_mpi.h"
#include <limits>
#include <vector>
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_doctest.h"
#include "pdm_dbbtree.h"

// #define HUGE_VAL 1.0e+30

MPI_TEST_CASE("[pdm_dbbtree] simple test",1) {

  double coord_tri_x[9] = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0};
  double coord_tri_y[9] = {0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0};
  double coord_tri_z[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int dim = 3;

  double tolerance = 1.e-3;

  const int         n_vtx            = 9;
  const int         n_cell           = 8;
  const int         n_tri_section_1  = 8;
  PDM_g_num_t connec_tri_1[24] = {1, 2, 5,
                                  1, 5, 4,
                                  2, 3, 6,
                                  2, 6, 5,
                                  4, 5, 8,
                                  4, 8, 7,
                                  5, 6, 9,
                                  5, 9, 8};

  std::vector<PDM_g_num_t> cell_ln_to_gn{1, 2, 3, 4, 5, 6, 7, 8};

  double* coords = (double *) malloc( 3 * n_vtx * sizeof(double));
  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {
    coords[3*i_vtx  ] = coord_tri_x[i_vtx];
    coords[3*i_vtx+1] = coord_tri_y[i_vtx];
    coords[3*i_vtx+2] = coord_tri_z[i_vtx];
  }

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  double global_extents[6] = {HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                             -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  double* extents = (double *) malloc( n_cell * 3 * 3 * sizeof(double));

  for(int i_cell = 0; i_cell < n_cell; ++i_cell) {

    // > Pour l'instant on prends la formule classique (see pdm_surf_mesh)
    for (int k1 = 0; k1 < 3; k1++) {
      extents[6*i_cell  +k1] =  std::numeric_limits<double>::max();
      extents[6*i_cell+3+k1] = -std::numeric_limits<double>::max();
    }

    for(int k = 3*i_cell; k < 3*(i_cell+1); ++k ){ // Car tri = Step by 3
      int i_vtx = connec_tri_1[k] - 1;
      for (int k1 = 0; k1 < 3; k1++) {
        extents[6*i_cell  +k1] = PDM_MIN (coords[3*i_vtx+k1], extents[6*i_cell  +k1]);
        extents[6*i_cell+3+k1] = PDM_MAX (coords[3*i_vtx+k1], extents[6*i_cell+3+k1]);
        // printf(" DEBUG k = %i | coord = %12.5e \n", k1, coords[3*i_vtx+k1]);
      }
    }

    double delta = -std::numeric_limits<double>::max();
    for (int k1 = 0; k1 < 3; k1++) {
      delta = PDM_MAX (delta, fabs (extents[6*i_cell+3+k1] - extents[6*i_cell  +k1]));
    }

    delta *= tolerance;

    for (int k1 = 0; k1 < 3; k1++) {
      extents[6*i_cell  +k1] +=  - delta;
      extents[6*i_cell+3+k1] +=    delta;
    }



    for (int k = 0; k < dim; k++) {
      global_extents[k]       = PDM_MIN(extents[2*dim*i_cell + k      ], global_extents[k]    );
      // global_extents[dim + k] = PDM_MAX(extents[(2*i_cell+1) * dim + k], global_extents[dim+k]);
      global_extents[dim + k] = PDM_MAX(extents[6*i_cell+3+k], global_extents[dim+k]);
    }
    // printf(" i_cell = %i | extent -> [", i_cell);
    // for (int k = 0; k < dim; k++) {
    //   printf("%12.5e %12.5e // ", extents[2*dim*i_cell + k], extents[6*i_cell+3+k]);
    // }
    // printf("] \n");

    // g_global_extents[]
  }

  // Print the final :
  // printf(" global_extents -> [");
  // for (int k = 0; k < dim; k++) {
  //   printf("%12.5e %12.5e // ", global_extents[k], global_extents[dim+k]);
  // }
  // printf("] \n");

  double max_range = -HUGE_VAL;
  double min_range = HUGE_VAL;

  for (int k = 0; k < dim; k++) {
    max_range = PDM_MAX(max_range, (global_extents[dim+k] - global_extents[k]));
    min_range = PDM_MIN(min_range, (global_extents[dim+k] - global_extents[k]));
  }

  for (int k = 0; k < dim; k++) {
    global_extents[k]     += -max_range * 1.1e-3; // On casse la symetrie !
    global_extents[dim+k] +=  max_range * 1e-3;
  }


  // Si parallèle il faut echangé le g_extents ... ( see pdm_overlay.c)
  PDM_dbbtree_t *dbbtreeA = PDM_dbbtree_create (pdm_comm, dim, global_extents);

  int n_part = 1;
  PDM_g_num_t** pcell_ln_to_gn = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t*));
  pcell_ln_to_gn[0] = cell_ln_to_gn.data();

  PDM_box_set_t  *boxesA = PDM_dbbtree_boxes_set (dbbtreeA,
                                                  n_part,
                                                  &n_cell,
                         (const double **)        &extents,
                         (const PDM_g_num_t **)   pcell_ln_to_gn);
  PDM_UNUSED(boxesA);
  PDM_UNUSED(n_tri_section_1);

  // Maintenant on souhaite localisé les points ...
  // int *boxesB_intersection_index;
  // int *boxesB_intersection_l_num;
  // PDM_box_set_t  *boxesB = PDM_dbbtree_intersect_boxes_set (dbbtreeA,
  //                                                           n_partB,
  //                                                           nEltsB,
  //                                                           extentsB,
  //                                                           gNumB,
  //                                                           &boxesB_intersection_index,
  //                                                           &boxesB_intersection_l_num);


  /*
   * Locate points
   */
  int         *pts_idx   = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  double      *pts_coord = NULL;

  int n_pts_pcloud = 1;
  int n_boxes      = n_cell;

  PDM_g_num_t pcloud_g_num[1] = {1};
  // double pcloud_coord[3] = {0.5, 0.5, 0.0001};
  double pcloud_coord[3] = {0.45, 0.45, 0.};


  PDM_dbbtree_points_inside_boxes (dbbtreeA,
                                   n_pts_pcloud,
                                   pcloud_g_num,
                                   pcloud_coord,
                                   n_boxes,
                                   cell_ln_to_gn.data(),
                                   &pts_idx,
                                   &pts_g_num,
                                   &pts_coord,
                                   0);


  // for(int i_box = 0; i_box < n_boxes; ++i_box) {
  //   printf("i_box = %i \n", i_box);
  //   for(int i = pts_idx[i_box]; i < pts_idx[i_box+1]; ++i) {
  //     printf(" | " PDM_FMT_G_NUM" %12.5e %12.5e %12.5e \n", pts_g_num[i], pts_coord[3*i], pts_coord[3*i+1], pts_coord[3*i+2]);
  //   }
  // }


  PDM_dbbtree_free (dbbtreeA);

  free(extents);
  free(coords);
}
