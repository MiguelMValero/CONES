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
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_closest_points.h"
#include "pdm_dcube_gen.h"
#include "pdm_geom_elem.h"
#include "pdm_version.h"
#include "pdm_logging.h"


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
     "  -c      <level>  Number of closest points (default : 10).\n\n"
     "  -t      <level>  Number of Target points (default : 10).\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_faceSeg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   nClosest Number of closest points
 * \param [inout]   nTgt     Number of Target points
 * \param [inout]   n_part   Number of partitions par process
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_faceSeg,
           double        *length,
           int           *nClosest,
           PDM_g_num_t   *nTgt,
           int           *n_part)
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
        long _n_faceSeg = atol(argv[i]);
        *n_faceSeg = (PDM_g_num_t) _n_faceSeg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *nClosest = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nTgt = atol(argv[i]);
        *nTgt = (PDM_g_num_t) _nTgt;
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
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



static void
_gen_cube_cell_centers
(
 PDM_MPI_Comm       comm,
 const PDM_g_num_t  n_faceSeg,
 const double       length,
 const double       zero_x,
 const double       zero_y,
 const double       zero_z,
 int               *npts,
 PDM_g_num_t      **g_num,
 double           **coord
 )
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *distribCell = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));

  PDM_g_num_t n_cell      = n_faceSeg * n_faceSeg * n_faceSeg;
  PDM_g_num_t n_face_face = n_faceSeg * n_faceSeg;

  // Define distribution
  distribCell[0] = 0;
  PDM_g_num_t stepCell = n_cell / n_rank;
  PDM_g_num_t remainderCell = n_cell % n_rank;

  for (int i = 1; i < n_rank + 1; i++) {
    distribCell[i]  = stepCell;
    const int i1 = i - 1;
    if (i1 < remainderCell)
      distribCell[i] += 1;
  }

  for (int i = 1; i < n_rank + 1; i++) {
    distribCell[i] += distribCell[i-1];
  }

  PDM_g_num_t _dn_cell = distribCell[i_rank+1] - distribCell[i_rank];

  const double step = length / (double) n_faceSeg;

  *g_num = malloc (sizeof(PDM_g_num_t) * _dn_cell);
  *coord = malloc (sizeof(double)      * _dn_cell * 3);

  int _npts = 0;
  for (PDM_g_num_t g = distribCell[i_rank]; g < distribCell[i_rank+1]; g++) {
    PDM_g_num_t i = g % n_faceSeg;
    PDM_g_num_t j = ((g - i) % n_face_face) / n_faceSeg;
    PDM_g_num_t k = (g - i - n_faceSeg * j) / n_face_face;

    (*coord)[3 * _npts    ] = (i + 0.5) * step + zero_x;
    (*coord)[3 * _npts + 1] = (j + 0.5) * step + zero_y;
    (*coord)[3 * _npts + 2] = (k + 0.5) * step + zero_z;
    (*g_num)[_npts++] = g + 1;//1 + i + n_faceSeg * j + n_face_face * k;
  }

  *npts = _npts;

  free (distribCell);
}





/**
 *
 * \brief  Main
 *
 */
// @@@param[n_proc] : 1,2,3,4
// @@@param[c] : 1,2,3,10
// @@@param[n] : 30, 60
// @@@param[t] : 10000, 20000
int main(int argc, char *argv[])
{
  int i_rank;
  int numProcs;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &numProcs);

  char *version = PDM_version_get();

  printf("Version de ParaDiGM : %s\n", version);
  free(version);

  /*
   *  Set default values
   */

  PDM_g_num_t  n_faceSeg = 10;
  double        length   = 1.;
  int           n_part   = 1;

  int n_closest_points = 10;
  PDM_g_num_t nTgt = 10;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_faceSeg,
             &length,
             &n_closest_points,
             &nTgt,
             &n_part);


  /* Define the target point cloud */
  srand(0);
  double      *tgt_coords = NULL;
  PDM_g_num_t *tgt_gnum   = NULL;
  int _n_tgt_l;

  double h = 0.5 * length / (double) n_faceSeg;
  PDM_point_cloud_gen_random (PDM_MPI_COMM_WORLD,
                              0, // seed
                              0, // geometric_g_num
                              nTgt,
                              h, h, h,
                              length - 2*h, length - 2*h, length - 2*h,
                              &_n_tgt_l,
                              &tgt_coords,
                              &tgt_gnum);

  /*
   *  Define the source point cloud (cell centers of cube)
   */
  n_part = 1;
  int _n_src_l;
  double *src_coords = NULL;
  PDM_g_num_t *src_gnum = NULL;
  _gen_cube_cell_centers (PDM_MPI_COMM_WORLD,
                          n_faceSeg,
                          length,
                          0.,
                          0.,
                          0.,
                          &_n_src_l,
                          &src_gnum,
                          &src_coords);

  n_closest_points = PDM_MIN (n_closest_points, n_faceSeg*n_faceSeg*n_faceSeg);


  /* Init closest points structure */
  PDM_closest_point_t* clsp = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                                         n_closest_points,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set (clsp,
                                       n_part,
                                       1);

  // set tgt point cloud
  PDM_closest_points_tgt_cloud_set (clsp,
                                    0,
                                    _n_tgt_l,
                                    tgt_coords,
                                    tgt_gnum);

  // set src point cloud
  PDM_closest_points_src_cloud_set (clsp,
                                    0,
                                    _n_src_l,
                                    src_coords,
                                    src_gnum);

  /* Compute closest points */
  PDM_closest_points_compute (clsp);

  PDM_closest_points_dump_times (clsp);

  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;

  PDM_closest_points_get (clsp,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);


  /* Check results */
  if (i_rank == 0) {
    printf("-- Check\n");
    fflush(stdout);
  }

  PDM_g_num_t n_wrong = 0;
  PDM_g_num_t n_tgt = 0;

  PDM_g_num_t *true_closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_closest_points);
  double      *true_closest_src_dist = malloc (sizeof(double) * n_closest_points);
  int n_cells_radius = (int) ceil(0.5 * pow((double) n_closest_points, 1./3.));

  int ijk0;
  int ijk_lo[3];
  int ijk_hi[3];
  double cell_side = length / ((double) n_faceSeg);
  double cell_ctr[3];

  for (int itgt = 0; itgt < _n_tgt_l; itgt++) {
    /*printf("[%d] %d/%d\n", i_rank, itgt, _n_tgt_l);*/

    n_tgt++;
    int wrong = 0;

    for (int l = 0; l < n_closest_points; l++) {
      true_closest_src_dist[l] = HUGE_VAL;
    }

    // find i,j,k of the cell that contains the target point
    // and define search region
    for (int idim = 0; idim < 3; idim++) {
      ijk0 = (int) floor(tgt_coords[3*itgt+idim] / cell_side);
      ijk0 = PDM_MIN (PDM_MAX (ijk0, 0), n_faceSeg-1);

      if (ijk0 < n_cells_radius) {
        ijk_lo[idim] = 0;
        ijk_hi[idim] = 2*n_cells_radius;
      } else if (ijk0 > n_faceSeg - n_cells_radius) {
        ijk_hi[idim] = n_faceSeg;
        ijk_lo[idim] = n_faceSeg - 2*n_cells_radius;
      } else {
        ijk_lo[idim] = ijk0 - n_cells_radius;
        ijk_hi[idim] = ijk0 + n_cells_radius;
      }
    }

    // inspect search region
    for (int k = ijk_lo[2]; k < ijk_hi[2]; k++) {
      cell_ctr[2] = (k + 0.5) * cell_side;
      for (int j = ijk_lo[1]; j < ijk_hi[1]; j++) {
        cell_ctr[1] = (j + 0.5) * cell_side;
        for (int i = ijk_lo[0]; i < ijk_hi[0]; i++) {
          cell_ctr[0] = (i + 0.5) * cell_side;

          PDM_g_num_t gnum = 1 + i + n_faceSeg*j + n_faceSeg*n_faceSeg*k;

          double dist = 0;
          for (int idim = 0; idim < 3; idim++) {
            double delta = tgt_coords[3*itgt+idim] - cell_ctr[idim];
            dist += delta * delta;
          }

          // insertion sort
          if (dist < true_closest_src_dist[n_closest_points-1]) {

            int l = n_closest_points - 1;
            while (l > 0 && dist < true_closest_src_dist[l-1]) {
              true_closest_src_gnum[l] = true_closest_src_gnum[l-1];
              true_closest_src_dist[l] = true_closest_src_dist[l-1];
              l--;
            }

            true_closest_src_gnum[l] = gnum;
            true_closest_src_dist[l] = dist;
          }
        }
      }
    }

    // check
    for (int l = 0; l < n_closest_points; l++) {
      /*printf("("PDM_FMT_G_NUM") %d/%d : "PDM_FMT_G_NUM"\t%f\n",
             tgt_gnum[itgt],
             l+1,
             n_closest_points,
             closest_src_gnum[n_closest_points*itgt + l],
             closest_src_dist[n_closest_points*itgt + l]);*/

      //if (closest_src_gnum[n_closest_points*itgt + l] != true_closest_src_gnum[l]) {
      /*if (closest_src_dist[n_closest_points*itgt + l] > true_closest_src_dist[l]) {
        printf("*** ERROR ("PDM_FMT_G_NUM") [%f %f %f]: %f / %f (relative err. = %f)\t\t"PDM_FMT_G_NUM" / "PDM_FMT_G_NUM"\n",
        tgt_gnum[itgt],
        tgt_coords[3*itgt], tgt_coords[3*itgt+1], tgt_coords[3*itgt+2],
        closest_src_dist[n_closest_points*itgt + l],
        true_closest_src_dist[l],
        (sqrt(closest_src_dist[n_closest_points*itgt + l]) - sqrt(true_closest_src_dist[l]))/sqrt(true_closest_src_dist[l]),
        closest_src_gnum[n_closest_points*itgt + l],
        true_closest_src_gnum[l]);
        }*/
      //assert (closest_src_gnum[n_closest_points*itgt + l] == true_closest_src_gnum[l]);
      //assert (closest_src_dist[n_closest_points*itgt + l] <= true_closest_src_dist[l]);

      if (closest_src_dist[n_closest_points*itgt + l] > true_closest_src_dist[l]) {
        printf("[%d] ("PDM_FMT_G_NUM") [%f %f %f]: %f / %f (relative err. = %f)\t\t"PDM_FMT_G_NUM" / "PDM_FMT_G_NUM"\n",
               i_rank,
               tgt_gnum[itgt],
               tgt_coords[3*itgt], tgt_coords[3*itgt+1], tgt_coords[3*itgt+2],
               closest_src_dist[n_closest_points*itgt + l],
               true_closest_src_dist[l],
               (sqrt(closest_src_dist[n_closest_points*itgt + l]) - sqrt(true_closest_src_dist[l]))/sqrt(true_closest_src_dist[l]),
               closest_src_gnum[n_closest_points*itgt + l],
               true_closest_src_gnum[l]);

        wrong = 1;
      }
    }

    if (wrong) {
      n_wrong++;
    }
  }
  free (true_closest_src_gnum);
  free (true_closest_src_dist);

  /*PDM_g_num_t wrong_percentage = 0;
  if (n_tgt > 0) {
    wrong_percentage = 100 * n_wrong / n_tgt;
  }
  printf("[%d] n_wrong = "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM" ("PDM_FMT_G_NUM"%%)\n",
         i_rank,
         n_wrong,
         n_tgt,
         wrong_percentage);*/

  PDM_g_num_t n_wrong_total, n_tgt_total;

  PDM_MPI_Reduce (&n_tgt,
                  &n_tgt_total,
                  1,
                  PDM__PDM_MPI_G_NUM,
                  PDM_MPI_SUM,
                  0,
                  PDM_MPI_COMM_WORLD);

  PDM_MPI_Reduce (&n_wrong,
                  &n_wrong_total,
                  1,
                  PDM__PDM_MPI_G_NUM,
                  PDM_MPI_SUM,
                  0,
                  PDM_MPI_COMM_WORLD);

  if (i_rank == 0) {
    PDM_g_num_t wrong_percentage_total = 0;
    if (n_tgt_total > 0) {
      wrong_percentage_total = 100 * n_wrong_total / n_tgt_total;
    }
    printf("n_wrong = "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM" ("PDM_FMT_G_NUM"%%)\n",
           n_wrong_total,
           n_tgt_total,
           wrong_percentage_total);

    //assert (n_wrong_total < 1);
  }







  /* Free */
  PDM_closest_points_free (clsp);

  free (tgt_coords);
  free (tgt_gnum);

  free (src_coords);
  free (src_gnum);


  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
