/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <sys/resource.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_timer.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"

#include "pdm_point_cloud_gen.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/


/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Generate a uniformly random point cloud inside a cuboid.
 *
 * \param [in]   comm                   MPI Communicator id
 * \param [in]   seed                   Random seed
 * \param [in]   gn_pts                 Global number of points in the cloud
 * \param [in]   geometric_g_num        Compute global ids from coordinates
 * \param [in]   x_min                  X-coordinate of the first cuboid corner
 * \param [in]   y_min                  Y-coordinate of the first cuboid corner
 * \param [in]   z_min                  Z-coordinate of the first cuboid corner
 * \param [in]   x_max                  X-coordinate of the opposite cuboid corner
 * \param [in]   y_max                  Y-coordinate of the opposite cuboid corner
 * \param [in]   z_max                  Z-coordinate of the opposite cuboid corner
 * \param [out]  ln_pts                 Local number of points in the cloud
 * \param [out]  coord                  XYZ-coordinates of the local points
 * \param [out]  g_num                  Global ids of the local points
 *
 */

void
PDM_point_cloud_gen_random
(
 PDM_MPI_Comm        comm,
 const int           seed,
 const int           geometric_g_num,
 const PDM_g_num_t   gn_pts,
 const double        x_min,
 const double        y_min,
 const double        z_min,
 const double        x_max,
 const double        y_max,
 const double        z_max,
 int                *ln_pts,
 double            **coord,
 PDM_g_num_t       **g_num
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t *distrib_pts = NULL;
  double      *dpts_coord = NULL;
  PDM_dpoint_cloud_gen_random(comm,
                              seed,
                              gn_pts,
                              x_min,
                              y_min,
                              z_min,
                              x_max,
                              y_max,
                              z_max,
                              &dpts_coord,
                              &distrib_pts);
  *coord = dpts_coord;

  *ln_pts = (int) (distrib_pts[i_rank+1] - distrib_pts[i_rank]);

  /**
   * Global numbers
   */
  if (geometric_g_num) {

    double length[3] = {x_max - x_min, y_max - y_min, z_max - z_min};

    double _char_length = 1e-6 * PDM_MAX(length[0], PDM_MAX(length[1], length[2]));

    double *char_length = malloc(sizeof(double) * (*ln_pts));

    for (int i = 0; i < *ln_pts; i++) {
      char_length[i] = _char_length;
    }

    PDM_gen_gnum_t* gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);

    PDM_gnum_set_from_coords (gen_gnum, 0, *ln_pts, *coord, char_length);

    PDM_gnum_compute (gen_gnum);

    *g_num = PDM_gnum_get (gen_gnum, 0);

    PDM_gnum_free (gen_gnum);
    free (char_length);
  }

  else {
    *g_num = malloc(sizeof(PDM_g_num_t) * (*ln_pts));
    for (int i = 0; i < *ln_pts; i++) {
      (*g_num)[i] = distrib_pts[i_rank] + i + 1;
    }
  }

  free (distrib_pts);
}

void
PDM_dpoint_cloud_gen_random
(
 PDM_MPI_Comm        comm,
 const int           seed,
 const PDM_g_num_t   gn_pts,
 const double        x_min,
 const double        y_min,
 const double        z_min,
 const double        x_max,
 const double        y_max,
 const double        z_max,
 double            **dpts_coord,
 PDM_g_num_t       **distrib_pts
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t *_distrib_pts = PDM_compute_uniform_entity_distribution(comm, gn_pts);

  /**
   * Coordinates
   */

  int dn_pts = (int) (_distrib_pts[i_rank+1] - _distrib_pts[i_rank]);
  double* _dpts_coord = malloc (dn_pts * 3 * sizeof(double));

  if (_dpts_coord == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Failed to allocate coords (size = %d * 3 * sizeof(double))\n", dn_pts);
  }

  double origin[3] = {x_min, y_min, z_min};

  double length[3] = {
    x_max - x_min,
    y_max - y_min,
    z_max - z_min
  };

  double i_rand_max = 1. / ((double) RAND_MAX);

  for (int i = 0; i < dn_pts; i++) {

    unsigned int _seed = (unsigned int) (_distrib_pts[i_rank] + i) + 1;
    srand(_seed + seed);

    for (int j = 0; j < 3; j++) {
      _dpts_coord[3*i + j] = origin[j] + length[j] * (double) rand() * i_rand_max;
    }
  }

  *dpts_coord  = _dpts_coord;
  *distrib_pts = _distrib_pts;

}




/**
 *
 * \brief Generate a cartesian point cloud inside a cuboid.
 *
 * \param [in]   comm                   MPI Communicator id
 * \param [in]   nx                     Number of points in X-direction
 * \param [in]   ny                     Number of points in Y-direction
 * \param [in]   nz                     Number of points in Z-direction
 * \param [in]   x_min                  X-coordinate of the first cuboid corner
 * \param [in]   y_min                  Y-coordinate of the first cuboid corner
 * \param [in]   z_min                  Z-coordinate of the first cuboid corner
 * \param [in]   x_max                  X-coordinate of the opposite cuboid corner
 * \param [in]   y_max                  Y-coordinate of the opposite cuboid corner
 * \param [in]   z_max                  Z-coordinate of the opposite cuboid corner
 * \param [out]  n_pts                  Local number of points in the cloud
 * \param [out]  pts_coord              XYZ-coordinates of the local points
 * \param [out]  pts_ln_to_gn           Global ids of the local points
 *
 */

void
PDM_point_cloud_gen_cartesian
(
 PDM_MPI_Comm        comm,
 const int           nx,
 const int           ny,
 const int           nz,
 const double        x_min,
 const double        y_min,
 const double        z_min,
 const double        x_max,
 const double        y_max,
 const double        z_max,
 int                *n_pts,
 double            **pts_coord,
 PDM_g_num_t       **pts_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t gn_pts = nx * ny * nz;

  PDM_g_num_t *distrib = PDM_compute_uniform_entity_distribution(comm,
                                                                 gn_pts);

  *n_pts = (int) (distrib[i_rank+1] - distrib[i_rank]);

  double      *_pts_coord    = malloc(sizeof(double)      * (*n_pts) * 3);
  PDM_g_num_t *_pts_ln_to_gn = malloc(sizeof(PDM_g_num_t) * (*n_pts));

  double step_x = (x_max - x_min) / (double) (nx - 1);
  double step_y = (y_max - y_min) / (double) (ny - 1);
  double step_z = (z_max - z_min) / (double) (nz - 1);

  for (int i = 0; i < *n_pts; ++i) {

    PDM_g_num_t g = distrib[i_rank] + i;

    _pts_ln_to_gn[i] = g + 1;

    PDM_g_num_t indi = g % nx;
    PDM_g_num_t indj = ((g - indi) / nx) % ny;
    PDM_g_num_t indk = g / (nx * ny);

    _pts_coord[3 * i    ] = x_min + indi * step_x;
    _pts_coord[3 * i + 1] = y_min + indj * step_y;
    _pts_coord[3 * i + 2] = z_min + indk * step_z;
  }

  free (distrib);


  *pts_coord    = _pts_coord;
  *pts_ln_to_gn = _pts_ln_to_gn;
}


void
PDM_dpoint_cloud_gen_cartesian
(
 PDM_MPI_Comm        comm,
 const int           nx,
 const int           ny,
 const int           nz,
 const double        x_min,
 const double        y_min,
 const double        z_min,
 const double        x_max,
 const double        y_max,
 const double        z_max,
 double            **dpts_coord,
 PDM_g_num_t       **distrib_pts
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t gn_pts = nx * ny * nz;

  PDM_g_num_t *distrib = PDM_compute_uniform_entity_distribution(comm,
                                                                 gn_pts);

  int dn_pts = (int) (distrib[i_rank+1] - distrib[i_rank]);

  double      *_pts_coord    = malloc(sizeof(double)      * dn_pts * 3);

  double step_x = 0;
  double step_y = 0;
  double step_z = 0;
  if (nx > 1) 
    step_x = (x_max - x_min) / (double) (nx - 1);
  if (ny > 1)
    step_y = (y_max - y_min) / (double) (ny - 1);
  if (nz > 1)
    step_z = (z_max - z_min) / (double) (nz - 1);

  for (int i = 0; i < dn_pts; ++i) {

    PDM_g_num_t g = distrib[i_rank] + i;

    PDM_g_num_t indi = g % nx;
    PDM_g_num_t indj = ((g - indi) / nx) % ny;
    PDM_g_num_t indk = g / (nx * ny);

    _pts_coord[3 * i    ] = x_min + indi * step_x;
    _pts_coord[3 * i + 1] = y_min + indj * step_y;
    _pts_coord[3 * i + 2] = z_min + indk * step_z;
  }

  *dpts_coord  = _pts_coord;
  *distrib_pts = distrib;
}
