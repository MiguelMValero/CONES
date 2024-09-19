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
#include "pdm_box_priv.h"
#include "pdm_box.h"
#include "pdm_box_tree.h"
#include "pdm_plane.h"
#include "pdm_linear_programming.h"

/*============================================================================
 * Static global variables
 *============================================================================*/

static const int    verbose = 0;
static const int    vtk     = 0;

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
  // return (double) rand() / (double) RAND_MAX;
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

/* Create a 4 plane volume */

static void
_create_volume_4planes
(
double  *edge,
double  *direction_pt,
double   theta,
double   eps,
double **n_in,
double **pt_plane_in
)
{
  *n_in        = malloc(sizeof(double) * 12);
  *pt_plane_in = malloc(sizeof(double) * 12);

  double *n = *n_in;
  double *pt_plane = *pt_plane_in;

  // Determine eps translation planes
  // B--{eps}--G
  double CB[3] = {edge[3]-edge[6], edge[4]-edge[7], edge[5]-edge[8]};
  pt_plane[0] = edge[3] + (1+eps) * CB[0];
  pt_plane[1] = edge[4] + (1+eps) * CB[1];
  pt_plane[2] = edge[5] + (1+eps) * CB[2];
  n[0] = -CB[0];
  n[1] = -CB[1];
  n[2] = -CB[2];

  // A--{eps}--H
  double CA[3] = {edge[3]-edge[0], edge[4]-edge[1], edge[5]-edge[2]};
  pt_plane[3] = edge[3] + (1+eps) * CA[0];
  pt_plane[4] = edge[4] + (1+eps) * CA[1];
  pt_plane[5] = edge[5] + (1+eps) * CA[2];
  n[3] = -CA[0];
  n[4] = -CA[1];
  n[5] = -CA[2];

  // Determine theta angle planes E---D---F
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);

  double AB[3] = {edge[6]-edge[0], edge[7]-edge[1], edge[8]-edge[2]};
  double inverse_module_AB = 1 / PDM_MODULE(AB);
  double AB_normalised[3] = {AB[0] * inverse_module_AB, AB[1] * inverse_module_AB, AB[2] * inverse_module_AB};

  pt_plane[6]  = (cos_theta + (1 - cos_theta) * AB_normalised[0] * AB_normalised[0]) * direction_pt[0];
  pt_plane[6] += (AB_normalised[0] * AB_normalised[1] * (1 - cos_theta) - AB_normalised[2] * sin_theta) * direction_pt[1];
  pt_plane[6] += (AB_normalised[0] * AB_normalised[2] * (1 - cos_theta) + AB_normalised[1] * sin_theta) * direction_pt[2];
  pt_plane[7]  = (AB_normalised[1] * AB_normalised[0] * (1 - cos_theta) + AB_normalised[2] * sin_theta) * direction_pt[0];
  pt_plane[7] += (cos_theta + (1 - cos_theta) * AB_normalised[1] * AB_normalised[1]) * direction_pt[1];
  pt_plane[7] += (AB_normalised[1] * AB_normalised[2] * (1 - cos_theta) - AB_normalised[0] * sin_theta) * direction_pt[2];
  pt_plane[8]  = (AB_normalised[2] * AB_normalised[0] * (1 - cos_theta) - AB_normalised[1] * sin_theta) * direction_pt[0];
  pt_plane[8] += (AB_normalised[2] * AB_normalised[1] * (1 - cos_theta) + AB_normalised[0] * sin_theta) * direction_pt[1];
  pt_plane[8] += (cos_theta + (1 - cos_theta) * AB_normalised[2] * AB_normalised[2]) * direction_pt[2];
  double prod_vect[3];
  double CE[3] = {pt_plane[6]-edge[3], pt_plane[7]-edge[4], pt_plane[8]-edge[5]};
  PDM_CROSS_PRODUCT(prod_vect, CA, CE);
  double ED[3] = {direction_pt[0] - pt_plane[6], direction_pt[1] - pt_plane[7], direction_pt[2] - pt_plane[8]};
  double sign = PDM_SIGN(PDM_DOT_PRODUCT(prod_vect, ED));
  n[6] = sign * prod_vect[0];
  n[7] = sign * prod_vect[1];
  n[8] = sign * prod_vect[2];

  double cos_minus_theta = cos(-theta);
  double sin_minus_theta = sin(-theta);

  pt_plane[9]   = (cos_minus_theta + (1 - cos_minus_theta) * AB_normalised[0] * AB_normalised[0]) * direction_pt[0];
  pt_plane[9]  += (AB_normalised[0] * AB_normalised[1] * (1 - cos_minus_theta) - AB_normalised[2] * sin_minus_theta) * direction_pt[1];
  pt_plane[9]  += (AB_normalised[0] * AB_normalised[2] * (1 - cos_minus_theta) + AB_normalised[1] * sin_minus_theta) * direction_pt[2];
  pt_plane[10]  = (AB_normalised[1] * AB_normalised[0] * (1 - cos_minus_theta) + AB_normalised[2] * sin_minus_theta) * direction_pt[0];
  pt_plane[10] += (cos_minus_theta + (1 - cos_minus_theta) * AB_normalised[1] * AB_normalised[1]) * direction_pt[1];
  pt_plane[10] += (AB_normalised[1] * AB_normalised[2] * (1 - cos_minus_theta) - AB_normalised[0] * sin_minus_theta) * direction_pt[2];
  pt_plane[11]  = (AB_normalised[2] * AB_normalised[0] * (1 - cos_minus_theta) - AB_normalised[1] * sin_minus_theta) * direction_pt[0];
  pt_plane[11] += (AB_normalised[2] * AB_normalised[1] * (1 - cos_minus_theta) + AB_normalised[0] * sin_minus_theta) * direction_pt[1];
  pt_plane[11] += (cos_minus_theta + (1 - cos_minus_theta) * AB_normalised[2] * AB_normalised[2]) * direction_pt[2];
  double CF[3] = {pt_plane[9]-edge[4], pt_plane[10]-edge[4], pt_plane[11]-edge[5]};
  PDM_CROSS_PRODUCT(prod_vect, CA, CF);
  double FD[3] = {direction_pt[0] - pt_plane[9], direction_pt[1] - pt_plane[10], direction_pt[2] -  pt_plane[11]};
  sign = PDM_SIGN(PDM_DOT_PRODUCT(prod_vect, FD));
  n[9]  = sign * prod_vect[0];
  n[10] = sign * prod_vect[1];
  n[11] = sign * prod_vect[2];
}

static void
_volume_4planes_to_4triangles
(
double       *edge,
double       *direction_pt,
double       *n,
double       *pt_plane,
double      **vtx_coord_in,
PDM_g_num_t **vtx_g_num_in,
int         **face_vtx_in
)
{
  *vtx_coord_in = malloc(sizeof(double) * 12 * 3);
  *vtx_g_num_in = malloc(sizeof(PDM_g_num_t) * 12);
  *face_vtx_in  = malloc(sizeof(int) * 12);

  double *vtx_coord = *vtx_coord_in;
  PDM_g_num_t *vtx_g_num = *vtx_g_num_in;
  int *face_vtx = *face_vtx_in;

  for (int i = 0; i < 12; i ++) {
    vtx_g_num[i] = i + 1;
    face_vtx[i]  = i +1;
  }

  double x[3];
  double origin[3];
  double normal[3];
  double proj_x[3];

  // translation plane 1
  // project E on 1
  for (int k = 0; k < 3; k++) {
    origin[k] = pt_plane[k]; // G
    normal[k] = n[k];
    x[k] = pt_plane[6 + k]; // E
    vtx_coord[k] = pt_plane[k]; // G
  }
  PDM_plane_projection(x, origin, normal, proj_x);
  // project D on 1
  for (int k = 0; k < 3; k++) {
    x[k] = direction_pt[k]; // D
    vtx_coord[3 + k] = proj_x[k]; // proj_E
  }
  PDM_plane_projection(x, origin, normal, proj_x);
  for (int k = 0; k < 3; k++) {
    vtx_coord[6 + k] = proj_x[k]; // proj_D
  }

  // translation plane 2
  // project E on 2
  for (int k = 0; k < 3; k++) {
    origin[k] = pt_plane[3 + k]; // H
    normal[k] = n[3 + k];
    x[k] = pt_plane[6 + k]; // E
    vtx_coord[9 + k] = pt_plane[3 + k]; // H
  }
  PDM_plane_projection(x, origin, normal, proj_x);
  // project D on 2
  for (int k = 0; k < 3; k++) {
    x[k] = direction_pt[k]; // D
    vtx_coord[12 + k] = proj_x[k]; // proj_E
  }
  PDM_plane_projection(x, origin, normal, proj_x);
  for (int k = 0; k < 3; k++) {
    vtx_coord[15 + k] = proj_x[k]; // proj_D
  }

  // rotation plane 3
  for (int k = 0; k < 3; k++) {
    vtx_coord[18 + 3*0 + k] = edge[k]; // A
    vtx_coord[18 + 3*1 + k] = edge[6 + k]; // B
    vtx_coord[18 + 3*2 + k] = pt_plane[6 + k]; // E
  }

  // rotation plane 4
  for (int k = 0; k < 3; k++) {
    vtx_coord[27 + 3*0 + k] = edge[k]; // A
    vtx_coord[27 + 3*1 + k] = edge[6 + k]; // B
    vtx_coord[27 + 3*2 + k] = pt_plane[9 + k]; // E
  }

}

/* box volume intersection function used in box_tree */

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

static int
_is_in_int_tab
(
int value,
int len_tab,
int *tab
)
{
  for (int i = 0; i < len_tab; i++) {
    if (tab[i] == value) {
      return 1;
    }
  }
  return 0;
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

  // Create volumes

  double *edge = malloc(sizeof(double) * 9) ; // A---C---B
  edge[0] = 0;
  edge[1] = 0;
  edge[2] = 0;
  edge[3] = 1;
  edge[4] = 0;
  edge[5] = 0;
  edge[6] = 2;
  edge[7] = 0;
  edge[8] = 0;
  double *direction_pt = malloc(sizeof(double) * 3) ; // C---D-->
  direction_pt[0] = 1;
  direction_pt[1] = 10;
  direction_pt[2] = 0;
  double theta = PDM_PI / 3;
  double eps   = 0;


  double *n        = NULL;
  double *pt_plane = NULL;

  _create_volume_4planes(edge, direction_pt, theta, eps, &n, &pt_plane);


  if (vtk) {
    double *vtx_coord =  NULL;
    PDM_g_num_t *vtx_g_num =  NULL;
    int *face_vtx = NULL;

    _volume_4planes_to_4triangles(edge,
                                  direction_pt,
                                  n,
                                  pt_plane,
                                  &vtx_coord,
                                  &vtx_g_num,
                                  &face_vtx);

    const char *filename = "planes.vtk";

    PDM_vtk_write_std_elements(filename,
                               12,
                               vtx_coord,
                               vtx_g_num,
                               PDM_MESH_NODAL_TRIA3,
                               4,
                               face_vtx,
                               NULL,
                               0,
                               NULL,
                               NULL);
    free(vtx_coord);
    free(vtx_g_num);
    free(face_vtx);
  }


  // Create boxes extents

  int gn_box = 10000;
  double      *box_extents  = NULL;
  PDM_g_num_t *box_ln_to_gn = NULL;

  int n_boxes = _generate_random_boxes(comm,
                                       gn_box,
                                       i_rank,
                                       &box_extents,
                                       &box_ln_to_gn);

  // Create boxes

  int *origin = PDM_array_zeros_int(3 * n_boxes);
  for (int i = 0 ; i < n_boxes; i++) {
    origin[3*i  ] = 0;
    origin[3*i+1] = 0;
    origin[3*i+2] = i;
  }

  PDM_box_set_t *boxes = PDM_box_set_create(3,
                                            1,
                                            0,
                                            n_boxes,
                                            box_ln_to_gn,
                                            box_extents,
                                            1, // n_part_orig
                                            &n_boxes, // n_boxes_orig
                                            origin, // origin
                                            comm);
  // Create PDM_box_tree_t
  PDM_box_tree_t* bt = PDM_box_tree_create(1000,
                                           10,
                                           10);

  // Build tree and associate boxes
  PDM_box_tree_set_boxes(bt,
                         boxes,
                         PDM_BOX_TREE_ASYNC_LEVEL);

  // WARNING: boxes outputed on internal nodes are wrong
  if (verbose) {
    PDM_box_tree_dump(bt);
  }

  // VTK boxes output
  if (vtk) {
    PDM_vtk_write_boxes("boxes.vtk",
                        n_boxes,
                        box_extents,
                        box_ln_to_gn);
  }

  // Brute force box intersection computation
  int *vol_boxes = malloc(sizeof(int) * n_boxes);
  int count = 0;
  int check;
  for (int j = 0; j < n_boxes; j++) {
    check = _box_intersect_volume(4, n, pt_plane, box_extents + 6*j);
    if (check) {
      vol_boxes[count++] = j;
    }
  }

  // Get for each volume the boxes it contains/intersects
  int *volume_box_idx   = NULL;
  int *volume_box_l_num = NULL;

  int n_volumes = 1;
  int *volume_plane_idx = malloc(sizeof(int) * (n_volumes +1));
  volume_plane_idx[0] = 0;
  volume_plane_idx[1] = 4;

  PDM_box_tree_intersect_volume_boxes(bt,
                                      -1,
                                      n_volumes,
                        (const int *) volume_plane_idx,
                                      n,
                                      pt_plane,
                                      &volume_box_idx,
                                      &volume_box_l_num);

  // Check if the same result as with brute force is obtained
  if (verbose) {
    log_trace("check if the result of brute force and the box tree is the same: \n");
    for (int i = 0; i < volume_box_idx[n_volumes]; i++) {
      if (!_is_in_int_tab(volume_box_l_num[i], count, vol_boxes)) {
        log_trace("box %d is not detected as intersecting by the box_tree\n", volume_box_l_num[i]);
      }
    }
  }

  free(vol_boxes);

  if (vtk) {
    PDM_box_tree_write_vtk("box_tree_normalized.vtk",
                           bt,
                           -1,
                           1);

    PDM_box_tree_write_vtk("box_tree.vtk",
                           bt,
                           -1,
                           0);
  }

  free(volume_box_idx);
  free(volume_box_l_num);
  free(volume_plane_idx);
  free(origin);
  free(box_ln_to_gn);
  free(box_extents);
  free(n);
  free(pt_plane);
  free(edge);
  free(direction_pt);
  PDM_box_set_destroy(&boxes);
  PDM_box_tree_destroy(&bt);

  PDM_MPI_Finalize ();
}
