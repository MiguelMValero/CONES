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
 char         **argv,
 int           *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else
      _usage (EXIT_FAILURE);
    i++;
  }
}

/*
 * \brief Determine if a box is on the side of the plane where the normal points to
 *
 * \param [in]  n                Normal vector (n = (a, b, c))
 * \param [in]  plane_pt         Point on plane (cartesian equation of plane being ax+by+cz+d=0 with d = -A.n)
 * \param [in]  box_extents      Box to determine
 *
 * \return 1 if the box is on the side of the plane where the normal points to, 0 otherwise
 */

static int
_plane_box_side
(
double *n,
double *plane_pt,
double *box_extents
)
{
  double box_pt[3];
  double vect[3];

  for (int x = 0; x < 4; x += 3) {
    box_pt[0] = box_extents[x];
    for (int y = 1; y < 5; y += 3) {
      box_pt[1] = box_extents[y];
      for (int z = 2; z < 6; z += 3) {
        box_pt[2] = box_extents[z];
        vect[0] = box_pt[0] - plane_pt[0]; vect[1] = box_pt[1] - plane_pt[1]; vect[2] = box_pt[2] - plane_pt[2];
        if (PDM_DOT_PRODUCT(vect, n) > 0) { // if >= 0 also considers when point is on the plane
          return 1;
        }
      } // end loop on z
    } // end loop on y
  } // end loop on x
  return 0;
}

/*
 * \brief Determine a box is in a given volume region
 *
 * \param [in]  n_planes      Number of planes difining the volume
 * \param [in]  n             Table of normal vector (n = (a, b, c)) of each considered plane
 * \param [in]  plane_pt      Table of a point on plane (cartesian equation of plane being ax+by+cz+d=0 with d = -A.n) for each considered plane
 * \param [in]  box_extents   Points to determine
 *
 * \return 1 if the plane_pt is on the side of the plane where the normal points to, 0 otherwise
 */


static int
_box_in_volume
(
int      n_planes,
double  *n,
double  *plane_pt,
double  *box_extents
)
{
  double n_iplane[3];
  double plane_pt_iplane[3];

  for (int iplane = 0; iplane < n_planes; iplane++) {

    n_iplane[0] = n[3*iplane];
    n_iplane[1] = n[3*iplane+1];
    n_iplane[2] = n[3*iplane+2];

    plane_pt_iplane[0] = plane_pt[3*iplane];
    plane_pt_iplane[1] = plane_pt[3*iplane+1];
    plane_pt_iplane[2] = plane_pt[3*iplane+2];

    if (_plane_box_side(n_iplane, plane_pt_iplane, box_extents) == 0) {
      return 0;
    }
  } // end loop on planes
  return 1;
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

  int visu = 0;
  _read_args (argc,
              argv,
              &visu);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &numProcs);

  /* Atomic test case */

  // Set up

  double edge[9] = {0, 0, 0, 1, 0, 0, 2, 0, 0}; // A---C---B
  double direction_pt[3] = {1, 20, 0}; // C---D-->
  double theta = PDM_PI / 3;
  double eps   = 1;

  double box_extents1[6] = {-1, 1, -1, 3, 5, 3}; // xmin, ymin, zmin, xmax, ymax, zmax
  double n[12];
  double pt_plane[12];

  // Determine eps translation planes
  // B--{eps}--G
  double CB[3] = {edge[3]-edge[6], edge[4]-edge[7], edge[5]-edge[8]};
  double inverse_module_CB = 1 / PDM_MODULE(CB);
  pt_plane[0] = edge[3] + (1+eps) * CB[0]; // * inverse_module_CB;
  pt_plane[1] = edge[4] + (1+eps) * CB[1]; // * inverse_module_CB;
  pt_plane[2] = edge[5] + (1+eps) * CB[2]; // * inverse_module_CB;
  n[0] = -CB[0];
  n[1] = -CB[1];
  n[2] = -CB[2];

  // A--{eps}--H
  double CA[3] = {edge[3]-edge[0], edge[4]-edge[1], edge[5]-edge[2]};
  double inverse_module_CA = 1 / PDM_MODULE(CA);
  pt_plane[3] = edge[3] + (1+eps) * CA[0]; // * inverse_module_CA;
  pt_plane[4] = edge[4] + (1+eps) * CA[1]; // * inverse_module_CA;
  pt_plane[5] = edge[5] + (1+eps) * CA[2]; // * inverse_module_CA;
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

  // Check if box is in volume
  int check = _box_in_volume(4, n, pt_plane, box_extents1);
  PDM_UNUSED(check);
  // log_trace("box is in volume = %d\n", check);

  // vtk output of atomic test case

  const char *filename1 = "box.vtk";
  PDM_g_num_t *box_g_num = malloc(sizeof(PDM_g_num_t) * 1);
  box_g_num[0] = 1;

  if (visu) {
    PDM_vtk_write_boxes(filename1,
                        1,
                        box_extents1,
                        box_g_num);
  }

  const char *filename2 = "line.vtk";
  double *coord = malloc(sizeof(double) * 6);
  coord[0] = edge[3];
  coord[1] = edge[4];
  coord[2] = edge[5];
  coord[3] = direction_pt[0];
  coord[4] = direction_pt[1];
  coord[5] = direction_pt[2];
  PDM_g_num_t *line_g_num = malloc(sizeof(PDM_g_num_t) * 1);
  line_g_num[0] = 1;

  if(visu) {
    PDM_vtk_write_lines(filename2,
                        1,
                        coord,
                        line_g_num,
                        NULL);
  }

  const char *filename3 = "planes.vtk";
  double *vtx_coord = malloc(sizeof(double) * 30);
  PDM_g_num_t *vtx_g_num = malloc(sizeof(PDM_g_num_t) * 10);
  int *face_vtx = malloc(sizeof(int) * 12);

  // A
  vtx_coord[0] = edge[0];
  vtx_coord[1] = edge[1];
  vtx_coord[2] = edge[2];

  // B
  vtx_coord[3] = edge[6];
  vtx_coord[4] = edge[7];
  vtx_coord[5] = edge[8];

  // E
  vtx_coord[6] = pt_plane[6];
  vtx_coord[7] = pt_plane[7];
  vtx_coord[8] = pt_plane[8];

  // F
  vtx_coord[9]  = pt_plane[9];
  vtx_coord[10] = pt_plane[10];
  vtx_coord[11] = pt_plane[11];

  // G
  vtx_coord[12] = pt_plane[0];
  vtx_coord[13] = pt_plane[1];
  vtx_coord[14] = pt_plane[2];

  // H
  vtx_coord[15] = pt_plane[3];
  vtx_coord[16] = pt_plane[4];
  vtx_coord[17] = pt_plane[5];

  double CH[3] = {pt_plane[3]-edge[3], pt_plane[4]-edge[4], pt_plane[5]-edge[5]};
  double module_CH = PDM_MODULE(CH);

  // I
  vtx_coord[18] = direction_pt[0] + module_CH * CA[0] * inverse_module_CA;
  vtx_coord[19] = direction_pt[1] + module_CH * CA[1] * inverse_module_CA;
  vtx_coord[20] = direction_pt[2] + module_CH * CA[2] * inverse_module_CA;

  // J
  vtx_coord[21] = pt_plane[6] + module_CH * CA[0] * inverse_module_CA;
  vtx_coord[22] = pt_plane[7] + module_CH * CA[1] * inverse_module_CA;
  vtx_coord[23] = pt_plane[8] + module_CH * CA[2] * inverse_module_CA;

  double CG[3] = {pt_plane[0]-edge[3], pt_plane[1]-edge[4], pt_plane[2]-edge[5]};
  double module_CG = PDM_MODULE(CG);

  // K
  vtx_coord[24] = direction_pt[0] + module_CG * CB[0] * inverse_module_CB;
  vtx_coord[25] = direction_pt[1] + module_CG * CB[1] * inverse_module_CB;
  vtx_coord[26] = direction_pt[2] + module_CG * CB[2] * inverse_module_CB;

  // L
  vtx_coord[27] = pt_plane[6] + module_CG * CB[0] * inverse_module_CB;
  vtx_coord[28] = pt_plane[7] + module_CG * CB[1] * inverse_module_CB;
  vtx_coord[29] = pt_plane[8] + module_CG * CB[2] * inverse_module_CB;

  for (int i = 0; i < 10; i++) {
    vtx_g_num[i] = i + 1;
  }

  face_vtx[0]  = 1;
  face_vtx[1]  = 2;
  face_vtx[2]  = 3;
  face_vtx[3]  = 1;
  face_vtx[4]  = 2;
  face_vtx[5]  = 4;
  face_vtx[6]  = 5;
  face_vtx[7]  = 9;
  face_vtx[8]  = 10;
  face_vtx[9]  = 6;
  face_vtx[10] = 7;
  face_vtx[11] = 8;

  if(visu) {
    PDM_vtk_write_std_elements(filename3,
                               10,
                               vtx_coord,
                               vtx_g_num,
                               PDM_MESH_NODAL_TRIA3,
                               4,
                               face_vtx,
                               NULL,
                               0,
                               NULL,
                               NULL);
  }

  const char *filename4 = "normal.vtk";

  double *vector_normal[1] = {n};

  const char* normal_name[] = {"n", 0};

  if(visu) {
    PDM_vtk_write_point_cloud_with_field(filename4,
                                         4,
                                         pt_plane,
                                         NULL,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL,
                                         1,
                         (const char **) &normal_name,
                       (const double **) &vector_normal,
                                         0,
                                         NULL,
                                         NULL);
  }

  /* Higher scale box-box test case */

  // double *n_box = malloc(sizeof(double) * 18);
  // double *pt_plane_box = malloc(sizeof(double) * 18);

  // _box_extents_to_plane(box_extents, n_box, pt_plane_box);

  // char *filename5 = "box_test.vtk";

  // double *vector_box[1] = {n_box};

  // const char* box_name[] = {"n", 0};

  // PDM_vtk_write_point_cloud_with_field(filename5,
  //                                      6,
  //                                      pt_plane_box,
  //                                      NULL,
  //                                      NULL,
  //                                      0,
  //                                      NULL,
  //                                      NULL,
  //                                      1,
  //                      (const char **) &box_name,
  //                    (const double **) &vector_box,
  //                                      0,
  //                                      NULL,
  //                                      NULL);

  // free(n_box);
  // free(pt_plane_box);

  PDM_g_num_t gn_box        = 10000;
  PDM_g_num_t gn_box_plane  = 10000;

  // Use random boxes generation ofpdm_t_intersect_line box

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

  //int  *box_tag[n_box_plane];
  double **box_tag = malloc(sizeof(double *) * n_box_plane);
  //char *box_tag_names[n_box_plane];
  char **box_tag_names = malloc(sizeof(char *) * n_box_plane);

  int check2;

  for (int i = 0; i < n_box_plane; i++) {

    _box_extents_iplane = box_extents_plane + 6 * i;

    _box_extents_to_plane(_box_extents_iplane, n_plane, pt_plane_box);

    for (int j = 0; j < n_box; j++) {

      _box_extents_current = box_extents + 6 * j;

      check = _box_in_volume(6, n_plane, pt_plane_box, _box_extents_current);
      check2 = _boxes_intersect_3d(_box_extents_iplane, _box_extents_current);
      if (check != check2) {
        log_trace("volume box %d and box box %d have check = %d but check2 = %d\n", i, j, check, check2);
      }
      box_intersects_box[j] = (double) check;
    }

    box_tag[i] = malloc(sizeof(double) * n_box);
    memcpy(box_tag[i], box_intersects_box, sizeof(double) * n_box);
    char tmp[90];
    sprintf(tmp, "intersects_box_%d", i);
    box_tag_names[i] = malloc(sizeof(char) * 90);
    strcpy(box_tag_names[i], tmp);


  }

  if(visu) {
    const char *filename6 = "scale_up_boxes.vtk";

    PDM_vtk_write_boxes_with_field(filename6,
                                   n_box,
                                   box_extents,
                                   NULL,
                                   n_box_plane,
                                   (const char **) box_tag_names,
                                   (const double **) box_tag);

    const char *filename7 = "scale_up_plane_boxes.vtk";

    PDM_vtk_write_boxes(filename7,
                        n_box_plane,
                        box_extents_plane,
                        NULL);
  }

  for (int i = 0; i < n_box; i++) {
    free(box_tag[i]);
    free(box_tag_names[i]);
  }
  free(coord);
  free(box_tag);
  free(box_tag_names);
  free(box_extents);
  free(box_ln_to_gn);
  free(box_extents_plane);
  free(box_ln_to_gn_plane);
  free(box_intersects_box);
  free(vtx_coord);
  free(vtx_g_num);
  free(face_vtx);
  free(line_g_num);
  free(box_g_num);
  PDM_MPI_Finalize ();

}
