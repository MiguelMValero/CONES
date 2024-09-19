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
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_linear_programming.h"

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const int    verbose = 0;
static const int    vtk     = 0;
static const int    only_lp = 1;

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
     "  -h               This message.\n\n");


  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else
      _usage (EXIT_FAILURE);
    i++;
  }
}

/*
 * \brief Determine planes from box extents
 *
 * \param [in]  box_extents   Input box caracteristics
 * \param [out]  n             Table of normal vector (n = (a, b, c)) of each considered plane
 * \param [out]  plane_pt      Table of a point on plane (cartesian equation of plane being ax+by+cz+d=0 with d = -A.n) for each considered plane
 *
 */

static void
_box_extents_to_plane
(
double *box_extents,
double *n,
double *pt_plane
)
{
  double vect1[3];
  double vect2[3];
  double vect3[3];
  double vectn[3];
  double sign;

  double xmin = box_extents[0];
  double ymin = box_extents[1];
  double zmin = box_extents[2];
  double xmax = box_extents[3];
  double ymax = box_extents[4];
  double zmax = box_extents[5];

  // normal x

  vect1[0] = xmax - xmin;
  vect1[1] = 0;
  vect1[2] = 0;

  vect2[0] = 0;
  vect2[1] = 0;
  vect2[2] = zmax - zmin;

  PDM_CROSS_PRODUCT(vectn, vect1, vect2);

  vect3[0] = 0;
  vect3[1] = ymax - ymin;
  vect3[2] = 0;

  sign = PDM_SIGN(PDM_DOT_PRODUCT(vect3, vectn));

  n[0] = sign * vectn[0];
  n[1] = sign * vectn[1];
  n[2] = sign * vectn[2];

  pt_plane[0] = xmin;
  pt_plane[1] = ymin;
  pt_plane[2] = zmin;

  n[3] = - sign * vectn[0];
  n[4] = - sign * vectn[1];
  n[5] = - sign * vectn[2];

  pt_plane[3] = xmin;
  pt_plane[4] = ymax;
  pt_plane[5] = zmin;

  // normal y

  vect1[0] = 0;
  vect1[1] = ymax - ymin;
  vect1[2] = 0;

  vect2[0] = 0;
  vect2[1] = 0;
  vect2[2] = zmax - zmin;

  PDM_CROSS_PRODUCT(vectn, vect1, vect2);

  vect3[0] = xmax - xmin;
  vect3[1] = 0;
  vect3[2] = 0;

  sign = PDM_SIGN(PDM_DOT_PRODUCT(vect3, vectn));

  n[6] = sign * vectn[0];
  n[7] = sign * vectn[1];
  n[8] = sign * vectn[2];

  pt_plane[6] = xmin;
  pt_plane[7] = ymin;
  pt_plane[8] = zmin;

  n[9]  = - sign * vectn[0];
  n[10] = - sign * vectn[1];
  n[11] = - sign * vectn[2];

  pt_plane[9]  = xmax;
  pt_plane[10] = ymin;
  pt_plane[11] = zmin;

  // normal z

  vect1[0] = xmax - xmin;
  vect1[1] = 0;
  vect1[2] = 0;

  vect2[0] = 0;
  vect2[1] = ymax - ymin;
  vect2[2] = 0;

  PDM_CROSS_PRODUCT(vectn, vect1, vect2);

  vect3[0] = 0;
  vect3[1] = 0;
  vect3[2] = zmax - zmin;

  sign = PDM_SIGN(PDM_DOT_PRODUCT(vect3, vectn));

  n[12] = sign * vectn[0];
  n[13] = sign * vectn[1];
  n[14] = sign * vectn[2];

  pt_plane[12] = xmin;
  pt_plane[13] = ymin;
  pt_plane[14] = zmin;

  n[15] = - sign * vectn[0];
  n[16] = - sign * vectn[1];
  n[17] = - sign * vectn[2];

  pt_plane[15] = xmin;
  pt_plane[16] = ymin;
  pt_plane[17] = zmax;

}

/* _random01 from pdm_t_intersect_line_box */

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

/* code from pdm_t_intersect_line_box */

static int
_generate_random_boxes
(
PDM_MPI_Comm  comm,
PDM_g_num_t   gn_box,
int           i_rank,
double      **box_extents,
PDM_g_num_t **box_ln_to_gn
)
{
  PDM_g_num_t *distrib_box = PDM_compute_uniform_entity_distribution (comm,
                                                                      gn_box);
  int n_box = (int) (distrib_box[i_rank+1] - distrib_box[i_rank]);
  for (PDM_g_num_t i = 0; i < 6*distrib_box[i_rank]; i++) {
    rand();
  }
  free (distrib_box);

  double *box_centers = malloc (sizeof(double) * n_box * 3);
  *box_extents = malloc (sizeof(double) * n_box * 6);
  double *_box_extents = *box_extents;
  for (int i = 0; i < n_box; i++) {
    for (int j = 0; j < 3; j++) {
      double x1 = _random01();
      double x2 = _random01();

      box_centers[3*i + j] = 0.5 * (x1 + x2);
      _box_extents[6*i + j]     = PDM_MIN (x1, x2);
      _box_extents[6*i + j + 3] = PDM_MAX (x1, x2);
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

  *box_ln_to_gn = PDM_gnum_get (gen_gnum, 0);

  PDM_gnum_free (gen_gnum);
  free (box_centers);

  return n_box;
}

/* _boxes_intersect_3d from pdm_box_tree */

inline static int
_boxes_intersect_3d
(
const double  *extents,
const double  *extentsB
)
{
  if (   extents[0] > extentsB[3] || extentsB[0] > extents[3]
   || extents[1] > extentsB[4] || extentsB[1] > extents[4]
   || extents[2] > extentsB[5] || extentsB[2] > extents[5])
    return 0;
  else
    return 1;
}

/*
 * \brief Determine if a box intersects a given volume region
 *
 * \param [in]  n_planes      Number of planes difining the volume
 * \param [in]  n             Table of normal vector (n = (a, b, c)) of each considered plane
 * \param [in]  plane_pt      Table of a point on plane (cartesian equation of plane being ax+by+cz+d=0 with d = -A.n) for each considered plane
 * \param [in]  box_extents   Points to determine
 *
 * \return 1 if the plane_pt is on the side of the plane where the normal points to, 0 otherwise
 */


static int
_box_intersect_volume
(
int      n_planes,
double  *n,
double  *plane_pt,
double  *box_extents
)
{
  double box_pt[3];
  double vect[3];
  double n_iplane[3];
  double plane_pt_iplane[3];

  int count_intersected_planes, count_points_not_intersect_plane;


  // All planes for one point
  for (int x = 0; x < 4; x += 3) {
    box_pt[0] = box_extents[x];
    for (int y = 1; y < 5; y += 3) {
      box_pt[1] = box_extents[y];
      for (int z = 2; z < 6; z += 3) {
        box_pt[2] = box_extents[z];

        count_intersected_planes = 0;

        for (int iplane = 0; iplane < n_planes; iplane++) {

          n_iplane[0] = n[3*iplane];
          n_iplane[1] = n[3*iplane+1];
          n_iplane[2] = n[3*iplane+2];

          plane_pt_iplane[0] = plane_pt[3*iplane];
          plane_pt_iplane[1] = plane_pt[3*iplane+1];
          plane_pt_iplane[2] = plane_pt[3*iplane+2];

          vect[0] = box_pt[0] - plane_pt_iplane[0]; vect[1] = box_pt[1] - plane_pt_iplane[1]; vect[2] = box_pt[2] - plane_pt_iplane[2];

          if (PDM_DOT_PRODUCT(vect, n_iplane) >= 0) {
            count_intersected_planes++;
          }

        } // end loop on planes

        if (count_intersected_planes == n_planes) {
          return 1;
        }

      }
    }
  }

  // All points for one plane
  for (int iplane = 0; iplane < n_planes; iplane++) {

    count_points_not_intersect_plane = 0;

    for (int x = 0; x < 4; x += 3) {
      box_pt[0] = box_extents[x];
      for (int y = 1; y < 5; y += 3) {
        box_pt[1] = box_extents[y];
        for (int z = 2; z < 6; z += 3) {
          box_pt[2] = box_extents[z];

          n_iplane[0] = n[3*iplane];
          n_iplane[1] = n[3*iplane+1];
          n_iplane[2] = n[3*iplane+2];

          plane_pt_iplane[0] = plane_pt[3*iplane];
          plane_pt_iplane[1] = plane_pt[3*iplane+1];
          plane_pt_iplane[2] = plane_pt[3*iplane+2];

          vect[0] = box_pt[0] - plane_pt_iplane[0]; vect[1] = box_pt[1] - plane_pt_iplane[1]; vect[2] = box_pt[2] - plane_pt_iplane[2];

          if (PDM_DOT_PRODUCT(vect, n_iplane) < 0) {
            count_points_not_intersect_plane++;
          }

        }
      }
    }

    if (count_points_not_intersect_plane == 6) {
      return 0;
    }

  }

  // Undefined case
  return PDM_lp_intersect_volume_box(n_planes, plane_pt, n, box_extents);
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init (&argc, &argv);

  int           i_rank;
  int           numProcs;

  _read_args (argc,
              argv);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &numProcs);

  PDM_g_num_t gn_box        = 1000;
  PDM_g_num_t gn_box_plane  = 1000;

  // Use random boxes generation of pdm_t_intersect_line box

  double *box_extents       = NULL;
  PDM_g_num_t *box_ln_to_gn = NULL;

  int n_box = _generate_random_boxes(comm,
                                     gn_box,
                                     i_rank,
                                     &box_extents,
                                     &box_ln_to_gn);

  double *box_extents_plane       = NULL;
  PDM_g_num_t *box_ln_to_gn_plane = NULL;

  int n_box_plane =_generate_random_boxes(comm,
                                          gn_box_plane,
                                          i_rank,
                                          &box_extents_plane,
                                          &box_ln_to_gn_plane);

  double *_box_extents_iplane = NULL;
  double *_box_extents_current = NULL;
  double n_plane[18];
  double pt_plane_box[18];

  double *box_intersects_box = malloc(sizeof(double) * n_box);

  double **box_tag = malloc(sizeof(double *) * n_box_plane);
  char **box_tag_names = malloc(sizeof(char *) * n_box_plane);

  int check1, check2;

  if (verbose) {
    log_trace("check if the result of classical box box intersection and lp is the same\n");
  }
  for (int i = 0; i < n_box_plane; i++) {

    _box_extents_iplane = box_extents_plane + 6 * i;

    _box_extents_to_plane(_box_extents_iplane, n_plane, pt_plane_box);

    for (int j = 0; j < n_box; j++) {

      _box_extents_current = box_extents + 6 * j;

      if (only_lp) {
        check1 = PDM_lp_intersect_volume_box(6, pt_plane_box, n_plane, _box_extents_current);
      } else {
        check1 = _box_intersect_volume(6, n_plane, pt_plane_box, _box_extents_current);
      }
      check2 = _boxes_intersect_3d(_box_extents_iplane, _box_extents_current);
      if (verbose) {
        if (check1 != check2) {
          log_trace("volume box %d and box box %d have check1 = %d but check2 = %d\n", i, j, check1, check2);
        }
      }
      box_intersects_box[j] = (double) check1;


    }

    box_tag[i] = malloc(sizeof(double) * n_box);
    memcpy(box_tag[i], box_intersects_box, sizeof(double) * n_box);
    char tmp[90];
    sprintf(tmp, "intersects_box_%d", i);
    box_tag_names[i] = malloc(sizeof(char) * 90);
    strcpy(box_tag_names[i], tmp);


  }

  if (vtk) {
    const char *filename6 = "boxes.vtk";

    PDM_vtk_write_boxes_with_field(filename6,
                                   n_box,
                                   box_extents,
                                   NULL,
                                   n_box_plane,
                                   (const char **) box_tag_names,
                                   (const double **) box_tag);

    const char *filename7 = "plane_boxes.vtk";

    PDM_vtk_write_boxes(filename7,
                        n_box_plane,
                        box_extents_plane,
                        NULL);
  }


  free(box_extents_plane);
  free(box_ln_to_gn_plane);
  free(box_extents);
  free(box_ln_to_gn);

  for (int j = 0; j < n_box; j++) {
    free(box_tag_names[j]);
    free(box_tag[j]);
  }
  free(box_tag_names);
  free(box_tag);
  free(box_intersects_box);

  PDM_MPI_Finalize ();

}
