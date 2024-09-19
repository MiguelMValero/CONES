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

#include "pdm_box_gen.h"

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
 * \brief Generate a random set of boxes
 *
 * \param [in]   comm                   MPI Communicator id
 * \param [in]   seed                   Random seed
 * \param [in]   geometric_g_num        Compute global ids from coordinates
 * \param [in]   gn_box                 Global number of boxes
 * \param [in]   min_size               Minimal box size
 * \param [in]   max_size               Maximal box size
 * \param [in]   x_min                  Minimal X-coordinate for box centers
 * \param [in]   y_min                  Minimal Y-coordinate for box centers
 * \param [in]   z_min                  Minimal Z-coordinate for box centers
 * \param [in]   x_max                  Maximal X-coordinate for box centers
 * \param [in]   y_max                  Maximal Y-coordinate for box centers
 * \param [in]   z_max                  Maximal Z-coordinate for box centers
 * \param [out]  n_box                  Local number of boxes
 * \param [out]  box_extents            Extents of the local boxes
 * \param [out]  box_ln_to_gn           Global ids of the local boxes
 *
 */

void
PDM_box_gen_random
(
 PDM_MPI_Comm   comm,
 int            seed,
 int            geometric_g_num,
 PDM_g_num_t    gn_box,
 double         min_size,
 double         max_size,
 double         x_min,
 double         y_min,
 double         z_min,
 double         x_max,
 double         y_max,
 double         z_max,
 int           *n_box,
 double       **box_extents,
 PDM_g_num_t  **box_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  double origin[3] = {x_min, y_min, z_min};
  double length[3] = {x_max - x_min, y_max - y_min, z_max - z_min};

  /*
   *  Generate random boxes
   */
  PDM_g_num_t *distrib_box = PDM_compute_uniform_entity_distribution (comm,
                                                                      gn_box);
  *n_box = (int) (distrib_box[i_rank+1] - distrib_box[i_rank]);

  *box_extents = malloc (sizeof(double) * (*n_box) * 6);
  for (int i = 0; i < (*n_box); i++) {

    unsigned int _seed = (unsigned int) (distrib_box[i_rank] + i) + 1;
    srand(_seed + seed);

    for (int j = 0; j < 3; j++) {
      double mid = origin[j] + length[j] * ((double) rand() / (double) RAND_MAX);
      double size = min_size + 0.5*(max_size - min_size) * ((double) rand() / (double) RAND_MAX);

      (*box_extents)[6*i + j]     = mid - size;
      (*box_extents)[6*i + j + 3] = mid + size;
    }
  }


  if (geometric_g_num) {
    double *box_centers = malloc (sizeof(double) * (*n_box) * 3);
    for (int i = 0; i < (*n_box); i++) {
      for (int j = 0; j < 3; j++) {
        box_centers[3*i+j] = 0.5 * ((*box_extents)[6*i+j] + (*box_extents)[6*i+j+3]);
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
                              *n_box,
                              box_centers,
                              NULL);

    PDM_gnum_compute (gen_gnum);
    free (box_centers);

    *box_ln_to_gn = PDM_gnum_get (gen_gnum, 0);

    PDM_gnum_free (gen_gnum);
  }
  else {
    *box_ln_to_gn = malloc(sizeof(PDM_g_num_t) * (*n_box));
    for (int i = 0; i < (*n_box); i++) {
      (*box_ln_to_gn)[i] = distrib_box[i_rank] + i + 1;
    }
  }
  free (distrib_box);

}


/**
 *
 * \brief Generate a cartesian set of boxes
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
 * \param [out]  n_box                  Local number of boxes
 * \param [out]  box_extents            Extents of the local boxes
 * \param [out]  box_ln_to_gn           Global ids of the local boxes
 *
 */

void
PDM_box_gen_cartesian
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
 int                *n_box,
 double            **box_extents,
 PDM_g_num_t       **box_ln_to_gn
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t gn_pts = (nx-1) * (ny-1) * (nz-1);

  PDM_g_num_t *distrib = PDM_compute_uniform_entity_distribution(comm,
                                                                 gn_pts);

  *n_box = (int) (distrib[i_rank+1] - distrib[i_rank]);

  double      *_box_extents  = malloc(sizeof(double)      * (*n_box) * 6);
  PDM_g_num_t *_box_ln_to_gn = malloc(sizeof(PDM_g_num_t) * (*n_box));

  double step_x = (x_max - x_min) / (double) (nx - 1);
  double step_y = (y_max - y_min) / (double) (ny - 1);
  double step_z = (z_max - z_min) / (double) (nz - 1);

  int n_box_x = nx - 1;
  int n_box_y = ny - 1;

  for (int i_box = 0; i_box < *n_box; ++i_box) {

    PDM_g_num_t g = distrib[i_rank] + i_box;

    _box_ln_to_gn[i_box] = g + 1;

    PDM_g_num_t ind_box_i = g % n_box_x;
    PDM_g_num_t ind_box_j = ((g - ind_box_i) / n_box_x) % n_box_y;
    PDM_g_num_t ind_box_k = g / (n_box_x * n_box_y);

    int i_vtx = ind_box_i + ind_box_j * nx + ind_box_k * nx * ny;

    PDM_g_num_t indi = i_vtx % nx;
    PDM_g_num_t indj = ((i_vtx - indi) / nx) % ny;
    PDM_g_num_t indk = i_vtx / (nx * ny);

    _box_extents[6 * i_box    ] = x_min + indi * step_x;
    _box_extents[6 * i_box + 1] = y_min + indj * step_y;
    _box_extents[6 * i_box + 2] = z_min + indk * step_z;

    _box_extents[6 * i_box + 3] = x_min + (indi+1) * step_x;
    _box_extents[6 * i_box + 4] = y_min + (indj+1) * step_y;
    _box_extents[6 * i_box + 5] = z_min + (indk+1) * step_z;

  }

  free (distrib);


  *box_extents  = _box_extents;
  *box_ln_to_gn = _box_ln_to_gn;
}
