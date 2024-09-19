#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_part_to_block.h"
#include "pdm_part_to_part.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"

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
     "  -t      <level>  Bounding boxes tolerance.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -p      <level>  Number of points to locate.\n\n"
     "  -octree          Use octree-based method.\n\n"
     "  -dbbree          Use dbbtree-based method.\n\n"
     "  -doctree          Use dbbtree-based method.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   tolerance  Bounding boxes tolerance
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_vtx_seg,
           double        *length,
           int           *rotation,
           double        *tolerance,
           double        *marge,
           int           *n_part,
           PDM_g_num_t   *n_pts,
           int           *post,
           int           *part_method,
           PDM_mesh_location_method_t *loc_method)
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
    else if (strcmp(argv[i], "-rot") == 0) {
        *rotation = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *marge = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_pts = atol(argv[i]);
        *n_pts = (PDM_g_num_t) _n_pts;
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_PART_SPLIT_HILBERT;
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = PDM_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-doctree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DOCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DBBTREE;
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void _rotate (const int  n_pts,
                     double    *coord)
{
  if (0) {
    double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
                      {0.3129918,  0.9447025, -0.0978434},
                      {-0.1593451,  0.1537920,  0.9751703}};

    for (int i = 0; i < n_pts; i++) {
      double x = coord[3*i];
      double y = coord[3*i+1];
      double z = coord[3*i+2];

      for (int j = 0; j < 3; j++) {
        coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
      }
    }
  }
  else {
    for (int i = 0; i < n_pts; i++) {
      double x = coord[3*i];
      coord[3*i+2] += 0.3 * cos(x * PDM_PI);
    }
  }
}

/**
 *
 * \brief  Random
 *
 */

static double
_frand_a_b(double a, double b){
  return ( rand()/(double)RAND_MAX ) * (b-a) + a;
}

static
void
_strip_cloud
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   n_vtx_seg,
 const double        xmin,
 const double        ymin,
 const double        zmin,
 const double        xmax,
 const double        ymax,
 const double        zmax,
 int                *ln_pts,
 double            **pts_coords_out,
 PDM_g_num_t       **g_num
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  double strip = 0.1;
  int n_pts   = ((n_vtx_seg * n_vtx_seg * n_vtx_seg) / n_rank);
  int n_pts_x = (int) (n_pts/(1.+(1.-4.*strip)+(1.-4.*strip) * (1.-4.*strip)));
  int n_pts_y = (int) ((1.-4.*strip) * n_pts);
  int n_pts_z = (int) ((1.-4.*strip) * (1.-4.*strip) * n_pts);

  if ((n_pts - ( n_pts_y + n_pts_z)) < 0)
    n_pts = n_pts_x + n_pts_y + n_pts_z;
  else
    n_pts_x = n_pts - ( n_pts_y + n_pts_z);

  *ln_pts = n_pts;
  *pts_coords_out  = malloc (sizeof(double) * 3 * n_pts);
  double *pts_coords  = *pts_coords_out;
  double *char_length = malloc (sizeof(double) * n_pts);

  double _char_length = (xmax-xmin)/n_vtx_seg;
  srand(i_rank+1); // Pas independant du parallèlisme mais pas de point double ...

  int idx  = 0;
  int idx2 = 0;
  for (int i = 0; i < n_pts_x/2; i++) {
    double x = _frand_a_b (xmin, xmin+strip);
    double y = _frand_a_b (ymin, ymax);
    double z = _frand_a_b (zmin, zmax);
    char_length[idx2++] = _char_length;
    pts_coords[idx++] = x;
    pts_coords[idx++] = y;
    pts_coords[idx++] = z;
  }

  for (int i =  n_pts_x/2; i < n_pts_x; i++) {
    double x = _frand_a_b (xmax-strip, xmax);
    double y = _frand_a_b (ymin, ymax);
    double z = _frand_a_b (zmin, zmax);
    char_length[idx2++] = _char_length;
    pts_coords[idx++] = x;
    pts_coords[idx++] = y;
    pts_coords[idx++] = z;
  }

  for (int i = 0; i < n_pts_y/2; i++) {
    double x = _frand_a_b (xmin+strip, xmax-strip);
    double y = _frand_a_b (ymin, ymin+strip);
    double z = _frand_a_b (zmin, zmax);
    char_length[idx2++] = _char_length;
    pts_coords[idx++] = x;
    pts_coords[idx++] = y;
    pts_coords[idx++] = z;
  }

  for (int i =  n_pts_y/2; i < n_pts_y; i++) {
    double x = _frand_a_b (xmin+strip, xmax-strip);
    double y = _frand_a_b (ymax-strip, ymax);
    double z = _frand_a_b (zmin, zmax);
    char_length[idx2++] = _char_length;
    pts_coords[idx++] = x;
    pts_coords[idx++] = y;
    pts_coords[idx++] = z;
  }

  for (int i = 0; i < n_pts_z/2; i++) {
    double x = _frand_a_b (xmin+strip, xmax-strip);
    double y = _frand_a_b (ymin+strip, ymax-strip);
    double z = _frand_a_b (zmin, zmin+strip);
    char_length[idx2++] = _char_length;
    pts_coords[idx++] = x;
    pts_coords[idx++] = y;
    pts_coords[idx++] = z;
  }

  for (int i =  n_pts_z/2; i < n_pts_z; i++) {
    double x = _frand_a_b (xmin+strip, xmax-strip);
    double y = _frand_a_b (ymin+strip, ymax-strip);
    double z = _frand_a_b (zmax-strip, zmax);
    char_length[idx2++] = _char_length;
    pts_coords[idx++] = x;
    pts_coords[idx++] = y;
    pts_coords[idx++] = z;
  }

  PDM_gen_gnum_t* gen_gnum_pts = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);

  PDM_gnum_set_from_coords (gen_gnum_pts, 0, n_pts, pts_coords, char_length);

  PDM_gnum_compute(gen_gnum_pts);

  *g_num =  PDM_gnum_get (gen_gnum_pts, 0);

  PDM_gnum_free(gen_gnum_pts);
  free(char_length);
}


static inline double
_eval_field
(
 double *xyz
 )
{
  return 1 + 2*xyz[0] + 3*xyz[1] + 4*xyz[2];
}



/**
 *
 * \brief  Main
 *
 */
// mpirun -np 2 ./test/pdm_t_mesh_location_dcube -n 3 -p 115 -doctree -post

// @@@param[n_proc] : 1,2,3,4
// @@@param[n] : 30,60
// @@@param[p] : 10000,20000
// @@@args[tree_kind] : -octree, -dbbree
int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t n_vtx_seg = 10;
  double      length    = 1.;
  int         rotation  = 0;
  double      tolerance = 1e-6;
  double      marge     = 0.;
  int         n_part    = 1;
  int         post      = 0;
  PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;

  PDM_g_num_t n_pts = 10;
  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &rotation,
             &tolerance,
             &marge,
             &n_part,
             &n_pts,
             &post,
             (int *) &part_method,
             &loc_method);


  /*
   *  Init
   */

  struct timeval t_elaps_debut;

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_rank      : %d\n", n_rank);
    PDM_printf ("  - n_vtx_seg   : "PDM_FMT_G_NUM"\n", n_vtx_seg);
    PDM_printf ("  - n_pts       : "PDM_FMT_G_NUM"\n", n_pts);
    PDM_printf ("  - length      : %f\n", length);
    PDM_printf ("  - tolerance   : %f\n", tolerance);
    PDM_printf ("  - part_method : %d\n", (int) part_method);
    PDM_printf ("  - loc_method  : %d\n", (int) loc_method);
  }

  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int          n_face_group;
  PDM_g_num_t *dface_cell = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx = NULL;
  double      *dvtx_coord = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group = NULL;
  int          dface_vtx_l;
  int          dface_group_l;

  /*
   *  Create distributed cube
   */

  const double xmin = 0;
  const double ymin = 0;
  const double zmin = 0;

  /*const double xmax = xmin + length;
    const double ymax = ymin + length;
    const double zmax = zmin + length;*/

  if (i_rank == 0) {
    printf("-- Build cube\n");
    fflush(stdout);
  }

  PDM_dcube_t* dcube = PDM_dcube_gen_init(PDM_MPI_COMM_WORLD,
                                          n_vtx_seg,
                                          length,
                                          xmin,
                                          ymin,
                                          zmin,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                        &n_face_group,
                        &dn_cell,
                        &dn_face,
                        &dn_vtx,
                        &dface_vtx_l,
                        &dface_group_l);

  PDM_dcube_gen_data_get(dcube,
                         &dface_cell,
                         &dface_vtx_idx,
                         &dface_vtx,
                         &dvtx_coord,
                         &dface_group_idx,
                         &dface_group);

  if (rotation) {
    _rotate (dn_vtx,
             dvtx_coord);
  }
  // int ppart_id = 0;

  gettimeofday(&t_elaps_debut, NULL);

  /*
   *  Create mesh partitions
   */

  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  if (i_rank == 0) {
    printf("-- Part\n");
    fflush(stdout);
  }

  PDM_part_t *ppart = PDM_part_create(PDM_MPI_COMM_WORLD,
                                      part_method,
                                      "PDM_PART_RENUM_CELL_NONE",
                                      "PDM_PART_RENUM_FACE_NONE",
                                      n_property_cell,
                                      renum_properties_cell,
                                      n_property_face,
                                      renum_properties_face,
                                      n_part,
                                      dn_cell,
                                      dn_face,
                                      dn_vtx,
                                      n_face_group,
                                      NULL,
                                      NULL,
                                      NULL,
                                      NULL,
                                      have_dcell_part,
                                      dcell_part,
                                      dface_cell,
                                      dface_vtx_idx,
                                      dface_vtx,
                                      NULL,
                                      dvtx_coord,
                                      NULL,
                                      dface_group_idx,
                                      dface_group);

  free(dcell_part);

  fflush(stdout);


  /************************
   *
   * Point cloud definition
   *
   ************************/
  if (i_rank == 0) {
    printf("-- Point cloud\n");
    fflush(stdout);
  }

  int n_pts_l;
  double      *pts_coords = NULL;
  PDM_g_num_t *pts_gnum   = NULL;

  marge *= length;
  int strip = 0;
  if(strip == 0) {
    // double llength = length/3;  // Desequilbrate
    double llength = length;

    PDM_point_cloud_gen_random (PDM_MPI_COMM_WORLD,
                                0, // seed
                                0, // geometric_g_num
                                n_pts,
                                -marge,
                                -marge,
                                -marge,
                                llength + marge,
                                llength + marge,
                                llength + marge,
                                &n_pts_l,
                                &pts_coords,
                                &pts_gnum);
  } else {

    const double xmax = xmin + length;
    const double ymax = ymin + length;
    const double zmax = zmin + length;
    _strip_cloud (PDM_MPI_COMM_WORLD,
                  n_vtx_seg,
                  xmin,
                  ymin,
                  zmin,
                  xmax,
                  ymax,
                  zmax,
                  &n_pts_l,
                  &pts_coords,
                  &pts_gnum);
  }

  if (rotation) {
    _rotate (n_pts_l,
             pts_coords);
  }

  if (post) {
    char filename[999];
    sprintf(filename, "point_cloud_%3.3d.vtk", i_rank);

    PDM_vtk_write_point_cloud (filename,
                               n_pts_l,
                               pts_coords,
                               pts_gnum,
                               NULL);
  }


  /************************
   *
   * Mesh location struct initializaiton
   *
   ************************/
  if (i_rank == 0) {
    printf("-- Mesh location_create\n");
    fflush(stdout);
  }

  PDM_mesh_location_t* mesh_loc = PDM_mesh_location_create (1,//const int n_point_cloud,
                                                            PDM_MPI_COMM_WORLD,
                                                            PDM_OWNERSHIP_KEEP);

  /* Set point cloud(s) */
  if (i_rank == 0) {
    printf("-- Mesh location_cloud_set\n");
    fflush(stdout);
  }
  PDM_mesh_location_n_part_cloud_set (mesh_loc,
                                      0,//i_point_cloud,
                                      1);//n_part

  PDM_mesh_location_cloud_set (mesh_loc,
                               0,//i_point_cloud,
                               0,//i_part,
                               n_pts_l,
                               pts_coords,
                               pts_gnum);

  if (i_rank == 0) {
    printf("-- Mesh location_mesh_set\n");
    fflush(stdout);
  }

  PDM_mesh_location_mesh_n_part_set (mesh_loc,
                                          n_part);

  /* Set mesh */
  for (int ipart = 0; ipart < n_part; ipart++) {

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_t_part;
    int s_cell_face;
    int s_face_vtx;
    int s_face_group;
    int n_edge_group2;

    PDM_part_part_dim_get(ppart,
                          ipart,
                          &n_cell,
                          &n_face,
                          &n_face_part_bound,
                          &n_vtx,
                          &n_proc,
                          &n_t_part,
                          &s_cell_face,
                          &s_face_vtx,
                          &s_face_group,
                          &n_edge_group2);

    int         *cell_tag;
    int         *cell_face_idx;
    int         *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int         *face_tag;
    int         *face_cell;
    int         *face_vtx_idx;
    int         *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int         *face_part_boundProcIdx;
    int         *face_part_boundPartIdx;
    int         *face_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *face_group_idx;
    int         *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart,
                           ipart,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_boundProcIdx,
                           &face_part_boundPartIdx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    PDM_mesh_location_part_set (mesh_loc,
                                ipart,
                                n_cell,
                                cell_face_idx,
                                cell_face,
                                cell_ln_to_gn,
                                n_face,
                                face_vtx_idx,
                                face_vtx,
                                face_ln_to_gn,
                                n_vtx,
                                vtx,
                                vtx_ln_to_gn);
  }

  /* Set location parameters */
  PDM_mesh_location_tolerance_set (mesh_loc,
                                   tolerance);

  PDM_mesh_location_method_set (mesh_loc,
                                loc_method);


  /* Compute location */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  // PDM_mesh_location_compute2 (mesh_loc);
  PDM_mesh_location_compute (mesh_loc);
  // PDM_mesh_location_compute(mesh_loc);

  PDM_mesh_location_dump_times (mesh_loc);

  int n_located = PDM_mesh_location_n_located_get (mesh_loc,
                                                   0,//i_point_cloud,
                                                   0);//i_part,

  int *located = PDM_mesh_location_located_get (mesh_loc,
                                                0,//i_point_cloud,
                                                0);

  int n_unlocated = PDM_mesh_location_n_unlocated_get (mesh_loc,
                                                       0,//i_point_cloud,
                                                       0);

  int *unlocated = PDM_mesh_location_unlocated_get (mesh_loc,
                                                    0,//i_point_cloud,
                                                    0);

  PDM_g_num_t *p_location    = NULL;
  double      *p_dist2  = NULL;
  double      *p_proj_coord  = NULL;
  PDM_mesh_location_point_location_get (mesh_loc,
                                        0,//i_point_cloud,
                                        0,//i_part,
                                        &p_location,
                                        &p_dist2,
                                        &p_proj_coord);



  if (1) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      int n_cell;
      int n_face;
      int n_face_part_bound;
      int n_vtx;
      int n_proc;
      int n_t_part;
      int s_cell_face;
      int s_face_vtx;
      int s_face_group;
      int n_edge_group2;

      PDM_part_part_dim_get(ppart,
                            ipart,
                            &n_cell,
                            &n_face,
                            &n_face_part_bound,
                            &n_vtx,
                            &n_proc,
                            &n_t_part,
                            &s_cell_face,
                            &s_face_vtx,
                            &s_face_group,
                            &n_edge_group2);

      int         *cell_tag;
      int         *cell_face_idx;
      int         *cell_face;
      PDM_g_num_t *cell_ln_to_gn;
      int         *face_tag;
      int         *face_cell;
      int         *face_vtx_idx;
      int         *face_vtx;
      PDM_g_num_t *face_ln_to_gn;
      int         *face_part_boundProcIdx;
      int         *face_part_boundPartIdx;
      int         *face_part_bound;
      int         *vtx_tag;
      double      *vtx;
      PDM_g_num_t *vtx_ln_to_gn;
      int         *face_group_idx;
      int         *face_group;
      PDM_g_num_t *face_group_ln_to_gn;

      PDM_part_part_val_get (ppart,
                             ipart,
                             &cell_tag,
                             &cell_face_idx,
                             &cell_face,
                             &cell_ln_to_gn,
                             &face_tag,
                             &face_cell,
                             &face_vtx_idx,
                             &face_vtx,
                             &face_ln_to_gn,
                             &face_part_boundProcIdx,
                             &face_part_boundPartIdx,
                             &face_part_bound,
                             &vtx_tag,
                             &vtx,
                             &vtx_ln_to_gn,
                             &face_group_idx,
                             &face_group,
                             &face_group_ln_to_gn);


      int *cell_vtx_idx;
      int *cell_vtx;
      PDM_mesh_location_cell_vertex_get (mesh_loc,
                                       ipart,//i_part_mesh,
                                       &cell_vtx_idx,
                                       &cell_vtx);

      int         *elt_pts_inside_idx;
      PDM_g_num_t *points_gnum;
      double      *points_coords;
      double      *points_uvw;
      int         *points_weights_idx;
      double      *points_weights;
      double      *points_dist2;
      double      *points_projected_coords;
      PDM_mesh_location_points_in_elt_get(mesh_loc,
                                          0,//i_point_cloud,
                                          ipart,
                                         &elt_pts_inside_idx,
                                         &points_gnum,
                                         &points_coords,
                                         &points_uvw,
                                         &points_weights_idx,
                                         &points_weights,
                                         &points_dist2,
                                         &points_projected_coords);


      if (0) {
        printf("cell_vtx : \n");
        for (int j = 0; j < n_cell; j++) {
          printf("%d :", j);
          for (int k = cell_vtx_idx[j]; k < cell_vtx_idx[j+1]; k++) {
            printf(" %d", cell_vtx[k]);
          }
          printf("\n");
        }
        printf("location in cell coords / uvw / dist2 / proj: \n");
        for (int j = 0; j < n_cell; j++) {
          printf(PDM_FMT_G_NUM" :", cell_ln_to_gn[j]);
          for (int k = elt_pts_inside_idx[j]; k < elt_pts_inside_idx[j+1]; k++) {
            printf(" "PDM_FMT_G_NUM" ooooo : %12.5e %12.5e %12.5e / %12.5e %12.5e %12.5e / %12.5e / %12.5e %12.5e %12.5e\n",
                   points_gnum[k],
                   points_coords[3*k],points_coords[3*k+1],points_coords[3*k+2],
                   points_uvw[3*k],points_uvw[3*k+1],points_uvw[3*k+2],points_dist2[k],
                   points_projected_coords[3*k],points_projected_coords[3*k+1],points_projected_coords[3*k+2]);
          // printf(" "PDM_FMT_G_NUM" ooooo : \n",
          //   points_gnum[k]);
          }
          printf("\n");
        }

        printf("location in cell weights : \n");
        for (int j = 0; j < n_cell; j++) {
          printf(PDM_FMT_G_NUM" :", cell_ln_to_gn[j]);
          for (int k = elt_pts_inside_idx[j]; k < elt_pts_inside_idx[j+1]; k++) {
            for (int k1 = points_weights_idx[k]; k1 < elt_pts_inside_idx[k+1]; k1++) {
              printf(" %12.5e", points_weights[k1]);
            }
            printf(" /");
          }
          printf("\n");
        }
      }
    }
  }

  // log_trace("n_pts_l = %d\n", n_pts_l);
  // PDM_log_trace_array_int(unlocated, n_unlocated, "unlocated : ");

  if (0) {
    printf("Unlocated %d :\n", n_unlocated);
    for (int k1 = 0; k1 < n_unlocated; k1++) {
      printf("%d\n", unlocated[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      printf("%d\n", located[k1]);
    }
    printf("\n");

    printf("Located %d :\n", n_located);
    for (int k1 = 0; k1 < n_located; k1++) {
      int ipt = located[k1] - 1;
      printf(PDM_FMT_G_NUM" : "PDM_FMT_G_NUM" / %12.5e %12.5e %12.5e / %12.5e / %12.5e %12.5e %12.5e",
        pts_gnum[ipt],  p_location[k1],
        pts_coords[3*ipt], pts_coords[3*ipt+1], pts_coords[3*ipt+2],
        p_dist2[k1],
        p_proj_coord[3*k1], p_proj_coord[3*k1+1], p_proj_coord[3*k1+2]);
      printf("\n");
    }
  }

  /* Check results */
  if (!rotation) {
    if (i_rank == 0) {
      printf("-- Check\n");
      fflush(stdout);
    }

    const double location_tolerance = 1.e-6 * length;

    const PDM_g_num_t n_cell_seg = n_vtx_seg - 1;
    const double cell_side = length / ((double) n_cell_seg);

    if(0 == 1) {
      printf("Unlocated :\n");
      for (int k1 = 0; k1 < n_unlocated; k1++) {
        printf("%d\n", unlocated[k1]);
      }
      printf("\n");
    }

    for (int k1 = 0; k1 < n_located; k1++) {
      int ipt = located[k1] - 1;
      double *p = pts_coords + 3*ipt;

      int i = (int) floor (p[0] / cell_side);
      int j = (int) floor (p[1] / cell_side);
      int k = (int) floor (p[2] / cell_side);

      PDM_g_num_t box_gnum = 1 + i + n_cell_seg*(j + n_cell_seg*k);

      if (p[0] < -tolerance || p[0] > length + tolerance ||
          p[1] < -tolerance || p[1] > length + tolerance ||
          p[2] < -tolerance || p[2] > length + tolerance) {
        box_gnum = -1;
      }

      //printf("%d: ("PDM_FMT_G_NUM") | ("PDM_FMT_G_NUM")\n", ipt, p_location[ipt], box_gnum);
      if (p_location[k1] != box_gnum) {
        // double *cp = p_proj_coord + 3*ipt;
        /*printf("%d ("PDM_FMT_G_NUM") (%.15lf %.15lf %.15lf): ("PDM_FMT_G_NUM") | ("PDM_FMT_G_NUM") proj : (%.15lf %.15lf %.15lf)\n",
               ipt, pts_gnum[ipt],
               p[0], p[1], p[2],
               p_location[ipt], box_gnum,
               cp[0], cp[1], cp[2]);
               printf("\n");*/

        //-->>
        double cell_min[3] = {cell_side * i,     cell_side * j,     cell_side * k};
        double cell_max[3] = {cell_side * (i+1), cell_side * (j+1), cell_side * (k+1)};
        /*printf("cell min = (%.15lf %.15lf %.15lf)\ncell max = (%.15lf %.15lf %.15lf)\n",
               cell_min[0], cell_min[1], cell_min[2],
               cell_max[0], cell_max[1], cell_max[2]);*/

        double dist = HUGE_VAL;
        for (int idim = 0; idim < 3; idim++) {
          double _dist1 = PDM_ABS (p[idim] - cell_min[idim]);
          double _dist2 = PDM_ABS (p[idim] - cell_max[idim]);
          double _dist = PDM_MIN (_dist1, _dist2);
          dist = PDM_MIN (dist, _dist);
        }
        //printf("distance = %e\n\n", dist);
        if (dist >= location_tolerance) {
          printf("point "PDM_FMT_G_NUM" : located in "PDM_FMT_G_NUM" (instead of "PDM_FMT_G_NUM"), dist = %e\n",
                 pts_gnum[ipt],
                 p_location[ipt],
                 box_gnum,
                 dist);
        }
        assert (dist < location_tolerance);
        //<<--
      }
    }
  }



  /*
   *  Check location (interpolation of an affine field)
   */
  double **src_field = malloc(sizeof(double *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_t_part;
    int s_cell_face;
    int s_face_vtx;
    int s_face_group;
    int n_edge_group2;

    PDM_part_part_dim_get(ppart,
                          ipart,
                          &n_cell,
                          &n_face,
                          &n_face_part_bound,
                          &n_vtx,
                          &n_proc,
                          &n_t_part,
                          &s_cell_face,
                          &s_face_vtx,
                          &s_face_group,
                          &n_edge_group2);

    int         *cell_tag;
    int         *cell_face_idx;
    int         *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int         *face_tag;
    int         *face_cell;
    int         *face_vtx_idx;
    int         *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int         *face_part_boundProcIdx;
    int         *face_part_boundPartIdx;
    int         *face_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *face_group_idx;
    int         *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart,
                           ipart,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_boundProcIdx,
                           &face_part_boundPartIdx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    src_field[ipart] = malloc(sizeof(double) * n_vtx);
    for (int i = 0; i < n_vtx; i++) {
      src_field[ipart][i] = _eval_field(&vtx[3*i]);
    }
  }

  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);

  double **send_field = malloc(sizeof(double *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    int         *elt_pts_idx        = NULL;
    PDM_g_num_t *elt_pts_gnum       = NULL;
    double      *elt_pts_coord      = NULL;
    double      *elt_pts_uvw        = NULL;
    int         *elt_pts_weight_idx = NULL;
    double      *elt_pts_weight     = NULL;
    double      *elt_pts_dist2      = NULL;
    double      *elt_pts_proj_coord = NULL;
    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        ipart,
                                        0, // i_point_cloud,
                                        &elt_pts_idx,
                                        &elt_pts_gnum,
                                        &elt_pts_coord,
                                        &elt_pts_uvw,
                                        &elt_pts_weight_idx,
                                        &elt_pts_weight,
                                        &elt_pts_dist2,
                                        &elt_pts_proj_coord);

    int *cell_vtx_idx = NULL;
    int *cell_vtx     = NULL;
    PDM_mesh_location_cell_vertex_get(mesh_loc,
                                      ipart,
                                      &cell_vtx_idx,
                                      &cell_vtx);

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_t_part;
    int s_cell_face;
    int s_face_vtx;
    int s_face_group;
    int n_edge_group2;

    PDM_part_part_dim_get(ppart,
                          ipart,
                          &n_cell,
                          &n_face,
                          &n_face_part_bound,
                          &n_vtx,
                          &n_proc,
                          &n_t_part,
                          &s_cell_face,
                          &s_face_vtx,
                          &s_face_group,
                          &n_edge_group2);

    send_field[ipart] = malloc(sizeof(double) * elt_pts_idx[n_cell]);
    for (int ielt = 0; ielt < n_cell; ielt++) {

      int *cv = cell_vtx + cell_vtx_idx[ielt];

      for (int idx_pt = elt_pts_idx[ielt]; idx_pt < elt_pts_idx[ielt+1]; idx_pt++) {
        send_field[ipart][idx_pt] = 0.;
        int idx_vtx = 0;
        for (int idx_w = elt_pts_weight_idx[idx_pt]; idx_w < elt_pts_weight_idx[idx_pt+1]; idx_w++) {
          int vtx_id = cv[idx_vtx++] - 1;
          send_field[ipart][idx_pt] += elt_pts_weight[idx_w] * src_field[ipart][vtx_id];
        }
      }

    }
  }


  double **recv_field = NULL;
  int request = -1;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) send_field,
                         NULL,
        (      void ***) &recv_field,
                         &request);

  PDM_part_to_part_iexch_wait(ptp, request);

  double lmax_err = 0.;
  for (int i = 0; i < n_located; i++) {
    int pt_id = located[i] - 1;

    double f = _eval_field(&pts_coords[3*pt_id]);

    double err = PDM_ABS(recv_field[0][i] - f);
    lmax_err = PDM_MAX(lmax_err, err);

    if (err > 1.e-12) {
      log_trace("point "PDM_FMT_G_NUM" (%f %f %f) located in elt "PDM_FMT_G_NUM" : error = %e (%20.16f / %20.16f)\n",
                pts_gnum[pt_id],
                pts_coords[3*pt_id], pts_coords[3*pt_id+1], pts_coords[3*pt_id+2],
                p_location[pt_id], err, recv_field[0][i], f);
    }
  }
  free(recv_field[0]);
  free(recv_field);

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(src_field [ipart]);
    free(send_field[ipart]);
  }
  free(src_field );
  free(send_field);


  double gmax_err;
  PDM_MPI_Allreduce(&lmax_err, &gmax_err, 1, PDM_MPI_DOUBLE,
                    PDM_MPI_MAX, PDM_MPI_COMM_WORLD);


  if (i_rank == 0) {
    printf("global max interpolation error = %e\n", gmax_err);
  }




  PDM_part_to_part_free(ptp);

  PDM_mesh_location_free (mesh_loc);

  PDM_part_free (ppart);

  free (pts_coords);
  free (pts_gnum);

  PDM_dcube_gen_free(dcube);
 
  double min_elaps_create;
  double max_elaps_create;
  double min_cpu_create;
  double max_cpu_create;
  double min_elaps_create2;
  double max_elaps_create2;
  double min_cpu_create2;
  double max_cpu_create2;
  double min_elaps_exch;
  double max_elaps_exch;
  double min_cpu_exch;
  double max_cpu_exch; 

  PDM_part_to_block_global_timer_get (PDM_MPI_COMM_WORLD,
                                      &min_elaps_create,
                                      &max_elaps_create,
                                      &min_cpu_create,
                                      &max_cpu_create,
                                      &min_elaps_create2,
                                      &max_elaps_create2,
                                      &min_cpu_create2,
                                      &max_cpu_create2,
                                      &min_elaps_exch,
                                      &max_elaps_exch,
                                      &min_cpu_exch,
                                      &max_cpu_exch);

  if (i_rank == 0) {
    printf("Global time in PDM_part_to_block : \n");
    printf("   - min max elaps create  : %12.5e %12.5e\n", min_elaps_create, max_elaps_create);
    printf("   - min max elaps create2 : %12.5e %12.5e\n", min_elaps_create2, max_elaps_create2);
    printf("   - min max elaps exch    : %12.5e %12.5e\n", min_elaps_exch, max_elaps_exch);
    fflush(stdout);
  }


  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
