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
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_distrib.h"
#include "pdm_dbbtree.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_binary_search.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"

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
     "  -b       <level> Number of boxes (default : 10).\n\n"
     "  -l       <level> Number of lines (default : 10).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nClosest   Number of closest points
 * \param [inout] nSrc   Number of Source points
 * \param [inout] nTgt   Number of Target points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *gn_box,
 PDM_g_num_t   *gn_line,
 int           *post
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-b") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long n = atol(argv[i]);
        *gn_box = (PDM_g_num_t) n;
      }
    }

    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long n = atol(argv[i]);
        *gn_line = (PDM_g_num_t) n;
      }
    }

    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}


/**
 *
 * \brief  Random value
 *
 *  \return a random double in [-1, 1]
 */

static double
_random01
(
 void
)
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / PDM_ABS (rsigna - rsignb);
  double resultat = sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}


inline static int
_intersect_line_box
(
 const int              dim,
 const double *restrict box_extents,
 const double *restrict line_origin,
 const double *restrict line_invdir
 )
{
  double tmin = 0.;
  double tmax = 1.;

  for (int i = 0; i < dim; i++) {
    double t1 = (box_extents[i]       - line_origin[i]) * line_invdir[i];
    double t2 = (box_extents[i + dim] - line_origin[i]) * line_invdir[i];

    if (line_invdir[i] < 0) {
      double tmp = t1;
      t1 = t2;
      t2 = tmp;
    }

    if (tmin > t2 || tmax < t1) {
      return 0;
    } else {
      tmin = PDM_MAX (tmin, t1);
      tmax = PDM_MIN (tmax, t2);
    }
  }

  return 1;
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

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t gn_box  = 10;
  PDM_g_num_t gn_line = 10;
  int         post    = 0;

  _read_args (argc,
              argv,
              &gn_box,
              &gn_line,
              &post);


  char filename[999];



  /*
   *  Generate random boxes
   */
  PDM_g_num_t *distrib_box = PDM_compute_uniform_entity_distribution (comm,
                                                                      gn_box);
  int n_box = (int) (distrib_box[i_rank+1] - distrib_box[i_rank]);
  for (PDM_g_num_t i = 0; i < 6*distrib_box[i_rank]; i++) {
    rand();
  }
  free (distrib_box);

  double *box_centers = malloc (sizeof(double) * n_box * 3);
  double *box_extents = malloc (sizeof(double) * n_box * 6);
  for (int i = 0; i < n_box; i++) {
    for (int j = 0; j < 3; j++) {
      double x1 = _random01();
      double x2 = _random01();

      box_centers[3*i + j] = 0.5 * (x1 + x2);
      box_extents[6*i + j]     = PDM_MIN (x1, x2);
      box_extents[6*i + j + 3] = PDM_MAX (x1, x2);
    }
  }


  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create (3,
                                              1,
                                              PDM_FALSE,
                                              1.e-3,
                                              comm,
                                              PDM_OWNERSHIP_USER);

  PDM_gnum_set_from_coords (gen_gnum,
                            0,
                            n_box,
                            box_centers,
                            NULL);

  PDM_gnum_compute (gen_gnum);

  PDM_g_num_t *box_ln_to_gn = PDM_gnum_get (gen_gnum, 0);

  PDM_gnum_free (gen_gnum);
  free (box_centers);




  /*
   *  Generate random lines
   */
  PDM_g_num_t *distrib_line = PDM_compute_uniform_entity_distribution (comm,
                                                                       gn_line);
  int n_line = (int) (distrib_line[i_rank+1] - distrib_line[i_rank]);
  for (PDM_g_num_t i = 0; i < 6*distrib_line[i_rank]; i++) {
    rand();
  }
  free (distrib_line);

  double *line_centers = malloc (sizeof(double) * n_line * 3);
  double *line_coords = malloc (sizeof(double) * n_line * 6);
  for (int i = 0; i < n_line; i++) {
    for (int j = 0; j < 3; j++) {
      double x1 = _random01();
      double x2 = _random01();

      line_centers[3*i + j] = 0.5 * (x1 + x2);
      line_coords[6*i + j]     = x1;
      line_coords[6*i + j + 3] = x2;
    }
  }


  gen_gnum = PDM_gnum_create (3,
                              1,
                              PDM_FALSE,
                              1.e-3,
                              comm,
                              PDM_OWNERSHIP_USER);

  PDM_gnum_set_from_coords (gen_gnum,
                            0,
                            n_line,
                            line_centers,
                            NULL);

  PDM_gnum_compute (gen_gnum);

  PDM_g_num_t *line_ln_to_gn = PDM_gnum_get (gen_gnum, 0);

  PDM_gnum_free (gen_gnum);
  free (line_centers);




  if (post) {
    sprintf(filename, "boxes_%2.2d.vtk", i_rank);
    PDM_vtk_write_boxes (filename,
                         n_box,
                         box_extents,
                         box_ln_to_gn);

    sprintf(filename, "lines_%2.2d.vtk", i_rank);
    PDM_vtk_write_lines (filename,
                         n_line,
                         line_coords,
                         line_ln_to_gn,
                         NULL);
  }




  /*
   *  Build dbbtree
   */
  const int dim = 3;
  double l_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL,
                         -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i = 0; i < n_box; i++) {
    for (int k = 0; k < 3; k++) {
      l_extents[k]     = PDM_MIN (l_extents[k],     box_extents[6*i + k]);
      l_extents[k + 3] = PDM_MAX (l_extents[k + 3], box_extents[6*i + k + 3]);
    }
  }

  double g_extents[6];
  PDM_MPI_Allreduce (l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce (l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX (max_range, g_extents[i+3] - g_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    g_extents[i]   -= max_range * 1.1e-3;
    g_extents[i+3] += max_range * 1.0e-3;
  }

  PDM_dbbtree_t *dbbt = PDM_dbbtree_create (comm, dim, g_extents);

  PDM_box_set_t *box_set = PDM_dbbtree_boxes_set (dbbt,
                                                  1,
                                                  &n_box,
                                                  (const double **) &box_extents,
                                                  (const PDM_g_num_t **) &box_ln_to_gn);


  /*
   *  Intersect lines with dbbtree
   */
  int         *intersecting_box_idx = NULL;
  PDM_g_num_t *intersecting_box_g_num = NULL;
  PDM_dbbtree_lines_intersect_boxes (dbbt,
                                     n_line,
                                     line_ln_to_gn,
                                     line_coords,
                                     &intersecting_box_idx,
                                     &intersecting_box_g_num);

  if (post) {
    for (int i = 0; i < n_line; i++) {
      log_trace("line "PDM_FMT_G_NUM": ", line_ln_to_gn[i]);
      for (int j = intersecting_box_idx[i]; j < intersecting_box_idx[i+1]; j++) {
        log_trace(PDM_FMT_G_NUM" ", intersecting_box_g_num[j]);
      }
      log_trace("\n");
    }
  }

  /*
   *  Check
   */
  PDM_g_num_t *all_box_ln_to_gn = malloc (sizeof(PDM_g_num_t) * gn_box);
  double      *all_box_extents  = malloc (sizeof(double)      * gn_box * 6);

  int *all_n_box = malloc (sizeof(int) * n_rank);
  PDM_MPI_Allgather (&n_box,    1, PDM_MPI_INT,
                     all_n_box, 1, PDM_MPI_INT,
                     comm);

  int *recv_shift = PDM_array_new_idx_from_sizes_int (all_n_box, n_rank);

  PDM_MPI_Allgatherv (box_ln_to_gn, n_box, PDM__PDM_MPI_G_NUM,
                      all_box_ln_to_gn, all_n_box, recv_shift, PDM__PDM_MPI_G_NUM,
                      comm);

  for (int i = 0; i < n_rank; i++) {
    all_n_box[i]    *= 6;
    recv_shift[i+1] *= 6;
  }

  PDM_MPI_Allgatherv (box_extents, 6*n_box, PDM_MPI_DOUBLE,
                      all_box_extents, all_n_box, recv_shift, PDM_MPI_DOUBLE,
                      comm);
  free (recv_shift);

  PDM_g_num_t ln_wrong = 0;
  double invdir[3];
  for (int iline = 0; iline < n_line; iline++) {

    double *orig = line_coords + 6*iline;
    double *dest = orig + 3;

    for (int i = 0; i < 3; i++) {
      double d = dest[i] - orig[i];
      if (PDM_ABS(d) < 1e-15) {
        invdir[i] = PDM_SIGN(d) * HUGE_VAL;
      } else {
        invdir[i] = 1. / d;
      }
    }

    int _n_intersect = intersecting_box_idx[iline+1] - intersecting_box_idx[iline];
    PDM_g_num_t *_intersect_g_num = intersecting_box_g_num + intersecting_box_idx[iline];

    for (PDM_g_num_t ibox = 0; ibox < gn_box; ibox++) {
      int intersect = _intersect_line_box (3,
                                           all_box_extents + 6*ibox,
                                           orig,
                                           invdir);

      int pos = -1;
      if (_n_intersect > 0) {
        pos = PDM_binary_search_long (all_box_ln_to_gn[ibox],
                                      _intersect_g_num,
                                      _n_intersect);
      }

      if ((pos >= 0) != intersect) {
        log_trace("error line "PDM_FMT_G_NUM", box "PDM_FMT_G_NUM" : pos = %d, intersect = %d\n", line_ln_to_gn[iline], all_box_ln_to_gn[ibox], pos, intersect);
        ln_wrong++;
      }
    }
  }
  free (all_n_box);
  free (all_box_ln_to_gn);
  free (all_box_extents);


  PDM_g_num_t gn_wrong;
  PDM_MPI_Allreduce (&ln_wrong, &gn_wrong, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  if (i_rank == 0) {
    printf(PDM_FMT_G_NUM" error(s)\n", gn_wrong);
  }


  /*
   *  Free memory
   */
  free (intersecting_box_idx);
  free (intersecting_box_g_num);

  PDM_dbbtree_free (dbbt);
  PDM_box_set_destroy (&box_set);

  free (box_extents);
  free (box_ln_to_gn);
  free (line_coords);
  free (line_ln_to_gn);

  PDM_MPI_Finalize();

  return gn_wrong;
}
