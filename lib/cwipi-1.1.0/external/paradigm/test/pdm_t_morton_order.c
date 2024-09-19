#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_morton.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_sort.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

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
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 int           *l,
 int           *lmax
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *l = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-lmax") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *lmax = atoi(argv[i]);
      }
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}



static uint64_t
_interleave
(
 const int               dimension,
 const PDM_morton_code_t c
 )
{
  int k = 0;
  uint64_t i = 0;
  for (PDM_morton_int_t l = 0; l < c.L; l++) {
    for (int j = dimension-1; j >= 0; j--) {
      uint64_t a = (c.X[j] >> l) & 1l;
      i += a << k;
      k++;
    }
  }

  return i;
}

static double
_code_to_double
(
 const int               dimension,
 const PDM_morton_code_t code
 )
{
  uint64_t i = _interleave(dimension, code);
  return (double) i;
}



inline static PDM_morton_code_t
_double_to_code(int     dim,
                double  input,
                int     level)
{
  int  l, child_id;
  PDM_morton_code_t  code;
  double coords[3] = {0.0, 0.0, 0.0};
  double l_mult = 1.0;

  const int max_level = 15; /* no more than 52 bits in mantissa / 3 */

  /* Build associated Morton code */

  code.L = max_level;

  if (input <= 0.0) {
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;
  }

  else if (input >= 1.0) {
    coords[0] = 1.0;
    coords[1] = 1.0;
    coords[2] = 1.0;
  }

  else if (dim == 3) {
    for (l = 0; l < max_level; l++) {
      l_mult *= 0.5;
      child_id = (int)(input*8);
      if (child_id > 7) child_id = 7;
      input = input*8 - child_id;
      coords[0] += child_id/4 * l_mult;
      coords[1] += (child_id%4)/2 * l_mult;
      coords[2] += child_id%2 * l_mult;
    }
  }

  else if (dim == 2) {
    coords[2] = 0;
    for (l = 0; l < max_level; l++) {
      l_mult *= 0.5;
      child_id = (int)(input*4);
      if (child_id > 3) child_id = 3;
      input = input*4 - child_id;
      coords[0] += child_id/2 * l_mult;
      coords[1] += child_id%2 * l_mult;
    }
  }

  else if (dim == 1) {
    coords[1] = 0;
    coords[2] = 0;
    for (l = 0; l < max_level; l++) {
      l_mult *= 0.5;
      child_id = (int)(input*2);
      if (child_id > 1) child_id = 1;
      input = input*2 - child_id;
      coords[0] += child_id * l_mult;
    }
  }

  code = PDM_morton_encode(dim, level, coords);

  return code;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
 int   argc,
 char *argv[]
 )
{
  PDM_MPI_Init (&argc, &argv);

  int level     = 3;
  int level_max = 10;
  _read_args(argc,
             argv,
             &level,
             &level_max);

  int n_seg = 1 << level;
  int n_pts = n_seg*n_seg*n_seg;

  double step = 1. / (double) n_seg;

  double *pts_coord = malloc(sizeof(double) * n_pts * 3);
  int idx = 0;
  // for (int k = 0; k < n_seg; k++) {
  //   for (int j = 0; j < n_seg; j++) {
  //     for (int i = 0; i < n_seg; i++) {
  //       pts_coord[idx++] = step*(i + 0.5);
  //       pts_coord[idx++] = step*(j + 0.5);
  //       pts_coord[idx++] = step*(k + 0.5);
  //     }
  //   }
  // }
  for (int i = 0; i < 3*n_pts; i++) {
    pts_coord[i] = (double) rand() / (double) RAND_MAX;
  }

  PDM_morton_code_t *pts_code = malloc(sizeof(PDM_morton_code_t) * n_pts);
  // for (int i = 0; i < n_pts; i++) {
  //   PDM_morton_code_t code = PDM_morton_encode(3,
  //                                              (PDM_morton_int_t) level,
  //                                              pts_coord + 3*i);

  //   PDM_morton_copy(code,
  //                   pts_code + i);
  // }
  double extents[6] = {0., 0., 0., 1., 1., 1.};
  double d[3], s[3];
  PDM_morton_encode_coords(3,
                           level_max,
                           extents,
                           n_pts,
                           pts_coord,
                           pts_code,
                           d,
                           s);

  PDM_morton_code_t *pts_code2 = malloc(sizeof(PDM_morton_code_t) * n_pts);
  memcpy(pts_code2, pts_code, sizeof(PDM_morton_code_t) * n_pts);

  PDM_morton_local_sort(n_pts,
                        pts_code2);


  int *order = malloc(sizeof(int) * n_pts);
  PDM_morton_local_order(n_pts,
                         pts_code,
                         order);


  // for (int i = 0; i < n_pts; i++) {
  //   log_trace("i = %d: L = %zu, X = %zu %zu %zu, order = %d\n",
  //             i,
  //             pts_code[i].L,
  //             pts_code[i].X[0], pts_code[i].X[1], pts_code[i].X[2],
  //             order[i]);
  // }

  int *reverse_order = malloc(sizeof(int) * n_pts);
  for (int i = 0; i < n_pts; i++) {
    reverse_order[order[i]] = i;
  }

  // uint64_t max = 1l << (level_max * 3);
  // log_trace("max = %zu\n", max);
  double normalization = (double) (1 << level_max);
  normalization = 1. / (normalization*normalization*normalization);
  // log_trace("normalization = %e\n", normalization);
  double *flat_code = malloc(sizeof(double) * n_pts);
  for (int i = 0; i < n_pts; i++) {
    flat_code[i] = _code_to_double(3, pts_code[i]);
    flat_code[i] *= normalization;
  }


  // PDM_vtk_write_point_cloud("morton_pts_coord.vtk",
  //                           n_pts,
  //                           pts_coord,
  //                           NULL,
  //                           reverse_order);



  // for (int i = 0; i < n_pts; i++) {
  //   int j = order[i];
  //   log_trace("i = %d: L = %zu, X = %zu %zu %zu / L = %zu, X = %zu %zu %zu, flat = %f\n",
  //             i,
  //             pts_code[j].L,
  //             pts_code[j].X[0], pts_code[j].X[1], pts_code[j].X[2],
  //             pts_code2[i].L,
  //             pts_code2[i].X[0], pts_code2[i].X[1], pts_code2[i].X[2],
  //             flat_code[j]);
  // }

  int *edge_vtx = malloc(sizeof(int) * 2 * (n_pts-1));
  for (int i = 0; i < n_pts-1; i++) {
    edge_vtx[2*i  ] = order[i  ] + 1;
    edge_vtx[2*i+1] = order[i+1] + 1;
  }

  free(pts_code2);

  double *_order = malloc(sizeof(double) * n_pts);
  for (int i = 0; i < n_pts; i++) {
    _order[i] = (double) reverse_order[i];
  }

  // const char   *field_name[]   = {"order", "flat_code"};
  // const double *field_value[2] = {_order, flat_code};



  // PDM_vtk_write_std_elements_ho_with_vtx_field("morton_curve.vtk",
  //                                              1,
  //                                              n_pts,
  //                                              pts_coord,
  //                                              NULL,
  //                                              PDM_MESH_NODAL_BAR2,
  //                                              n_pts-1,
  //                                              edge_vtx,
  //                                              NULL,
  //                                              0,
  //                                              NULL,
  //                                              NULL,
  //                                              2,
  //                                              field_name,
  //                                              field_value);

  free(edge_vtx);


  double *box_extents = malloc(sizeof(double) * n_pts * 6);
  idx = 0;
  for (int k = 0; k < n_seg; k++) {
    for (int j = 0; j < n_seg; j++) {
      for (int i = 0; i < n_seg; i++) {
        box_extents[idx++] = i*step;
        box_extents[idx++] = j*step;
        box_extents[idx++] = k*step;
        box_extents[idx++] = (i+1)*step;
        box_extents[idx++] = (j+1)*step;
        box_extents[idx++] = (k+1)*step;
      }
    }
  }

  PDM_g_num_t *g_num = malloc(sizeof(PDM_g_num_t) * n_pts);
  for (int i = 0; i < n_pts; i++) {
    g_num[i] = (PDM_g_num_t) reverse_order[i];
  }

  // PDM_vtk_write_boxes("morton_grid.vtk",
  //                     n_pts,
  //                     box_extents,
  //                     g_num);
  free(box_extents);
  free(g_num);


  int *order2 = malloc(sizeof(int) * n_pts);
  for (int i = 0; i < n_pts; i++) {
    order2[i] = i;
  }
  PDM_sort_double(flat_code,
                  order2,
                  n_pts);


  // for (int i = 0; i < n_pts; i++) {
  //   log_trace("%6d / %6d\n", order[i], order2[i]);
  // }
  free(order2);


  for (int i = 0; i < n_pts; i++) {
    PDM_morton_code_t code = _double_to_code(3, flat_code[i], level_max);
    double            dble = _code_to_double(3, code);
    dble *= normalization;
    // log_trace("i = %d: L = %zu, X = %zu %zu %zu / L = %zu, X = %zu %zu %zu, flat = %f / %f\n",
    //           i,
    //           pts_code[i].L,
    //           pts_code[i].X[0], pts_code[i].X[1], pts_code[i].X[2],
    //           code.L,
    //           code.X[0], code.X[1], code.X[2],
    //           flat_code[i], dble);
  }

  free(pts_coord);
  free(pts_code);
  free(order);
  free(reverse_order);
  free(_order);
  free(flat_code);

  PDM_MPI_Finalize ();

  return 0;

}
