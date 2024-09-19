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
#include "pdm_mesh_location.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_part_to_part.h"
#include "pdm_vtk.h"
#include "pdm_writer.h"

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
_read_args
(
 int                         argc,
 char                      **argv,
 PDM_g_num_t                *n_mesh,
 PDM_g_num_t                *n_cloud,
 double                     *radius,
 double                     *tolerance,
 int                        *n_part_mesh,
 int                        *n_part_cloud,
 int                        *post,
 int                        *part_method,
 PDM_mesh_location_method_t *loc_method
 )
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
        long _n = atol(argv[i]);
        *n_mesh  = (PDM_g_num_t) _n;
        *n_cloud = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nm") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n = atol(argv[i]);
        *n_mesh = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-np") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n = atol(argv[i]);
        *n_cloud = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-r") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *radius = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part_mesh  = atoi(argv[i]);
        *n_part_cloud = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part_mesh") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part_mesh = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_part_cloud") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part_cloud = atoi(argv[i]);
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


static void
_rotate
(
 const int     n_pts,
       double *coord
 )
{
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


static void
_dump_mesh
(
 const char     *name,
 PDM_MPI_Comm    comm,
 int             n_part,
 int            *pn_vtx,
 double        **pvtx_coord,
 PDM_g_num_t   **pvtx_ln_to_gn,
 int            *pn_face,
 int           **pface_vtx_idx,
 int           **pface_vtx,
 PDM_g_num_t   **pface_ln_to_gn,
 const int       n_elt_field,
 const char    **elt_field_name,
 double       ***elt_field_value,
 const int       n_vtx_field,
 const char    **vtx_field_name,
 double       ***vtx_field_value
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        name,
                                        name,
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       name,
                                       n_part);

  int id_var_num_part = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var_elt[n_elt_field];
  for (int i = 0; i < n_elt_field; i++) {
    id_var_elt[i] = PDM_writer_var_create(wrt,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          elt_field_name[i]);
  }

  int id_var_vtx[n_vtx_field];
  for (int i = 0; i < n_vtx_field; i++) {
    id_var_vtx[i] = PDM_writer_var_create(wrt,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_VERTICES,
                                          vtx_field_name[i]);
  }


  PDM_writer_step_beg(wrt, 0.);

  int id_block = PDM_writer_geom_bloc_add(wrt,
                                          id_geom,
                                          PDM_WRITER_POLY_2D,
                                          PDM_OWNERSHIP_USER);

  PDM_real_t **val_num_part = malloc(sizeof(PDM_real_t  *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_writer_geom_coord_set(wrt,
                             id_geom,
                             ipart,
                             pn_vtx[ipart],
                             pvtx_coord[ipart],
                             pvtx_ln_to_gn[ipart],
                             PDM_OWNERSHIP_USER);

    PDM_writer_geom_bloc_poly2d_set(wrt,
                                    id_geom,
                                    id_block,
                                    ipart,
                                    pn_face[ipart],
                                    pface_vtx_idx[ipart],
                                    pface_vtx[ipart],
                                    pface_ln_to_gn[ipart]);

    val_num_part[ipart] = malloc(sizeof(PDM_real_t) * pn_face[ipart]);
    for (int i = 0; i < pn_face[ipart]; i++) {
      val_num_part[ipart][i] = n_part*i_rank + ipart;
    }

    PDM_writer_var_set(wrt,
                       id_var_num_part,
                       id_geom,
                       ipart,
                       val_num_part[ipart]);

    for (int i = 0; i < n_elt_field; i++) {
      PDM_writer_var_set(wrt,
                         id_var_elt[i],
                         id_geom,
                         ipart,
                         elt_field_value[i][ipart]);
    }

    for (int i = 0; i < n_vtx_field; i++) {
      PDM_writer_var_set(wrt,
                         id_var_vtx[i],
                         id_geom,
                         ipart,
                         vtx_field_value[i][ipart]);
    }
  }

  PDM_writer_geom_write(wrt,
                        id_geom);

  PDM_writer_var_write(wrt,
                       id_var_num_part);
  PDM_writer_var_free(wrt,
                      id_var_num_part);

  for (int i = 0; i < n_elt_field; i++) {
    PDM_writer_var_write(wrt,
                         id_var_elt[i]);
    PDM_writer_var_free(wrt,
                        id_var_elt[i]);
  }

  for (int i = 0; i < n_vtx_field; i++) {
    PDM_writer_var_write(wrt,
                         id_var_vtx[i]);
    PDM_writer_var_free(wrt,
                        id_var_vtx[i]);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(val_num_part[ipart]);
  }
  free(val_num_part);

  PDM_writer_step_end(wrt);

  PDM_writer_free(wrt);
}




static void
_dump_point_cloud
(
 const char     *name,
 PDM_MPI_Comm    comm,
 int             n_part,
 int            *pn_pts,
 double        **ppts_coord,
 PDM_g_num_t   **ppts_ln_to_gn,
 const int       n_field,
 const char    **field_name,
 double       ***field_value
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        name,
                                        name,
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       name,
                                       n_part);

  int id_var_num_part = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var[n_field];
  for (int i = 0; i < n_field; i++) {
    id_var[i] = PDM_writer_var_create(wrt,
                                      PDM_WRITER_OFF,
                                      PDM_WRITER_VAR_SCALAR,
                                      PDM_WRITER_VAR_ELEMENTS,
                                      field_name[i]);
  }

  PDM_writer_step_beg(wrt, 0.);

  int id_block = PDM_writer_geom_bloc_add(wrt,
                                          id_geom,
                                          PDM_WRITER_POINT,
                                          PDM_OWNERSHIP_USER);

  PDM_real_t **val_num_part = malloc(sizeof(PDM_real_t  *) * n_part);
  int **connec = malloc(sizeof(int *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              ipart,
                              pn_pts[ipart],
                              ppts_coord[ipart],
                              ppts_ln_to_gn[ipart],
                              PDM_OWNERSHIP_USER);

    connec[ipart] = malloc(sizeof(int) * pn_pts[ipart]);
    for (int i = 0; i < pn_pts[ipart]; i++) {
      connec[ipart][i] = i+1;
    }
    PDM_writer_geom_bloc_std_set(wrt,
                                 id_geom,
                                 id_block,
                                 ipart,
                                 pn_pts[ipart],
                                 connec[ipart],
                                 ppts_ln_to_gn[ipart]);

    val_num_part[ipart] = malloc(sizeof(PDM_real_t) * pn_pts[ipart]);
    for (int i = 0; i < pn_pts[ipart]; i++) {
      val_num_part[ipart][i] = n_part*i_rank + ipart;
    }

    PDM_writer_var_set(wrt,
                       id_var_num_part,
                       id_geom,
                       ipart,
                       val_num_part[ipart]);

    for (int i = 0; i < n_field; i++) {
      PDM_writer_var_set(wrt,
                         id_var[i],
                         id_geom,
                         ipart,
                         field_value[i][ipart]);
    }
  }

  PDM_writer_geom_write(wrt,
                        id_geom);

  PDM_writer_var_write(wrt,
                       id_var_num_part);
  PDM_writer_var_free(wrt,
                      id_var_num_part);

  for (int i = 0; i < n_field; i++) {
    PDM_writer_var_write(wrt,
                         id_var[i]);
    PDM_writer_var_free(wrt,
                        id_var[i]);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {
    free(val_num_part[ipart]);
    free(connec[ipart]);
  }
  free(val_num_part);
  free(connec);

  PDM_writer_step_end(wrt);

  PDM_writer_free(wrt);
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t n_mesh       = 5;
  PDM_g_num_t n_cloud      = 10;
  double      radius       = 1.;
  double      tolerance    = 1e-3;
  int         n_part_mesh  = 1;
  int         n_part_cloud = 1;
  int         post         = 0;
  PDM_split_dual_t     part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_mesh,
              &n_cloud,
              &radius,
              &tolerance,
              &n_part_mesh,
              &n_part_cloud,
              &post,
      (int *) &part_method,
              &loc_method);

  tolerance *= radius;

  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Create partitionned surface mesh
   */

  if (i_rank == 0) {
    printf("-- Generate surface mesh\n");
    fflush(stdout);
  }

  int          *pn_vtx         = NULL;
  double      **pvtx_coord     = NULL;
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;
  int          *pn_face        = NULL;
  int         **pface_vtx_idx  = NULL;
  int         **pface_vtx      = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;
  PDM_sphere_surf_icosphere_gen_part(comm,
                                     n_mesh,
                                     0.,
                                     0.,
                                     0.,
                                     radius,
                                     n_part_mesh,
                                     part_method,
                                     &pn_vtx,
                                     &pvtx_coord,
                                     &pvtx_ln_to_gn,
                                     &pn_face,
                                     &pface_vtx_idx,
                                     &pface_vtx,
                                     &pface_ln_to_gn);


  /*
   *  Create partitionned point cloud
   */

  if (i_rank == 0) {
    printf("-- Generate point cloud\n");
    fflush(stdout);
  }

  int          *pn_pts          = NULL;
  double      **ppts_coord      = NULL;
  PDM_g_num_t **ppts_ln_to_gn   = NULL;
  int          *pn_face2        = NULL;
  int         **pface_vtx_idx2  = NULL;
  int         **pface_vtx2      = NULL;
  PDM_g_num_t **pface_ln_to_gn2 = NULL;
  PDM_sphere_surf_icosphere_gen_part(comm,
                                     n_cloud,
                                     0.,
                                     0.,
                                     0.,
                                     radius,
                                     n_part_cloud,
                                     part_method,
                                     &pn_pts,
                                     &ppts_coord,
                                     &ppts_ln_to_gn,
                                     &pn_face2,
                                     &pface_vtx_idx2,
                                     &pface_vtx2,
                                     &pface_ln_to_gn2);

  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    _rotate(pn_pts[ipart],
            ppts_coord[ipart]);
  }



  /*
   *  Mesh location
   */
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create(1,
                                                           comm,
                                                           PDM_OWNERSHIP_KEEP);

  /* Set point cloud */
  PDM_mesh_location_n_part_cloud_set(mesh_loc,
                                     0,//i_point_cloud,
                                     n_part_cloud);

  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    PDM_mesh_location_cloud_set(mesh_loc,
                                0,//i_point_cloud,
                                ipart,
                                pn_pts[ipart],
                                ppts_coord[ipart],
                                ppts_ln_to_gn[ipart]);
  }

  /* Set mesh */
  PDM_part_mesh_nodal_t *mesh_nodal = PDM_part_mesh_nodal_create(2,
                                                                 n_part_mesh,
                                                                 comm);

  int i_section = PDM_part_mesh_nodal_section_add(mesh_nodal,
                                                  PDM_MESH_NODAL_TRIA3);//POLY_2D);

  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    PDM_part_mesh_nodal_coord_set(mesh_nodal,
                                  ipart,
                                  pn_vtx[ipart],
                                  pvtx_coord[ipart],
                                  pvtx_ln_to_gn[ipart],
                                  PDM_OWNERSHIP_USER);

    int         *parent_num          = NULL;
    PDM_g_num_t *parent_entity_g_num = NULL;
    PDM_part_mesh_nodal_section_std_set(mesh_nodal,
                                        ipart,
                                        i_section,
                                        pn_face[ipart],
                                        pface_vtx[ipart],
                                        pface_ln_to_gn[ipart],
                                        parent_num,
                                        parent_entity_g_num,
                                        PDM_OWNERSHIP_USER);
  }

  PDM_mesh_location_shared_nodal_mesh_set(mesh_loc,
                                          mesh_nodal);


  /* Set location parameters */
  PDM_mesh_location_tolerance_set(mesh_loc,
                                  tolerance);

  PDM_mesh_location_method_set(mesh_loc,
                               loc_method);

  /*
   * Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }


  // PDM_mesh_location_compute(mesh_loc);
  PDM_mesh_location_compute(mesh_loc);




  /*
   *  Exchange data between point cloud and mesh
   */
  int                 **pelt_pts_idx        = malloc(sizeof(int         *) * n_part_mesh);
  PDM_g_num_t         **pelt_pts_gnum       = malloc(sizeof(PDM_g_num_t *) * n_part_mesh);
  double              **pelt_pts_coord      = malloc(sizeof(double      *) * n_part_mesh);
  double              **pelt_pts_uvw        = malloc(sizeof(double      *) * n_part_mesh);
  int                 **pelt_pts_weight_idx = malloc(sizeof(int         *) * n_part_mesh);
  double              **pelt_pts_weight     = malloc(sizeof(double      *) * n_part_mesh);
  double              **pelt_pts_dist2      = malloc(sizeof(double      *) * n_part_mesh);
  double              **pelt_pts_proj_coord = malloc(sizeof(double      *) * n_part_mesh);

  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        ipart,
                                        0, // i_point_cloud,
                                        &pelt_pts_idx       [ipart],
                                        &pelt_pts_gnum      [ipart],
                                        &pelt_pts_coord     [ipart],
                                        &pelt_pts_uvw       [ipart],
                                        &pelt_pts_weight_idx[ipart],
                                        &pelt_pts_weight    [ipart],
                                        &pelt_pts_dist2     [ipart],
                                        &pelt_pts_proj_coord[ipart]);
  }

  /* Create exchange protocol */
  // PDM_part_to_part_t *ptp = PDM_part_to_part_create((const PDM_g_num_t **) pface_ln_to_gn,
  //                                                   (const int          *) pn_face,
  //                                                   n_part_mesh,
  //                                                   (const PDM_g_num_t **) ppts_ln_to_gn,
  //                                                   (const int          *) pn_pts,
  //                                                   n_part_cloud,
  //                                                   (const int         **) pelt_pts_idx,
  //                                                   (const PDM_g_num_t **) pelt_pts_gnum,
  //                                                   comm);
  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);

  int  *n_located;
  int **located;
  PDM_part_to_part_ref_lnum2_get(ptp,
                                 &n_located,
                                 &located);

  int  *n_unlocated;
  int **unlocated;
  PDM_part_to_part_unref_lnum2_get(ptp,
                                   &n_unlocated,
                                   &unlocated);

  int         **gnum1_come_from_idx;
  PDM_g_num_t **gnum1_come_from;
  PDM_part_to_part_gnum1_come_from_get(ptp,
                                       &gnum1_come_from_idx,
                                       &gnum1_come_from);

  int          *n_elt1;
  int         **part1_to_part2_idx;
  PDM_g_num_t **part1_to_part2;
  PDM_part_to_part_part1_to_part2_get(ptp,
                                      &n_elt1,
                                      &part1_to_part2_idx,
                                      &part1_to_part2);

  // for (int ipart = 0; ipart < n_part_cloud; ipart++) {
  //   int _n_located = PDM_mesh_location_n_located_get(mesh_loc,
  //                                                    0,
  //                                                    ipart);

  //   int _n_unlocated = PDM_mesh_location_n_unlocated_get(mesh_loc,
  //                                                        0,
  //                                                        ipart);

  //   int *_located = PDM_mesh_location_located_get(mesh_loc,
  //                                                 0,
  //                                                 ipart);

  //   int *_unlocated = PDM_mesh_location_unlocated_get(mesh_loc,
  //                                                     0,
  //                                                     ipart);

  //   PDM_log_trace_array_int(located[ipart], n_located[ipart], "located (ptp) : ");
  //   PDM_log_trace_array_int(_located,       _n_located,       "located (ml)  : ");

  //   PDM_log_trace_array_int(unlocated[ipart], n_unlocated[ipart], "unlocated (ptp) : ");
  //   PDM_log_trace_array_int(_unlocated,       _n_unlocated,       "unlocated (ml)  : ");

  //   log_trace("n_elt1 = %d / %d\n", n_elt1[ipart], pn_face[ipart]);
  //   PDM_log_trace_array_int(part1_to_part2_idx[ipart], n_elt1 [ipart], "elt_pts_idx (ptp) : ");
  //   PDM_log_trace_array_int(pelt_pts_idx      [ipart], pn_face[ipart], "elt_pts_idx (ml)  : ");
  // }



  /*
   *  Exchange mesh (cell-centered) -> point cloud
   */
  double **pface_field1 = malloc(sizeof(double *) * n_part_mesh);
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    pface_field1[ipart] = malloc(sizeof(double) * pn_face[ipart]);

    for (int i = 0; i < pn_face[ipart]; i++) {
      double y = 0.;

      for (int idx = pface_vtx_idx[ipart][i]; idx < pface_vtx_idx[ipart][i+1]; idx++) {
        int vtx_id = pface_vtx[ipart][idx] - 1;
        y += pvtx_coord[ipart][3*vtx_id+1];
      }

      if (pface_vtx_idx[ipart][i+1] > pface_vtx_idx[ipart][i]) {
        y /= (double) (pface_vtx_idx[ipart][i+1] - pface_vtx_idx[ipart][i]);
      }

      pface_field1[ipart][i] = sin(5*y);
    }
  }

  double **pface_field1_p1p2 = malloc(sizeof(double *) * n_part_mesh);
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    pface_field1_p1p2[ipart] = malloc(sizeof(double) * pelt_pts_idx[ipart][pn_face[ipart]]);

    for (int i = 0; i < pn_face[ipart]; i++) {
      for (int idx = pelt_pts_idx[ipart][i]; idx < pelt_pts_idx[ipart][i+1]; idx++) {
        pface_field1_p1p2[ipart][idx] = pface_field1[ipart][i];
      }
    }
  }

  int request_elt_pts;
  double **tmp_ppts_field = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) pface_field1,
        // (const void  **) pface_field1_p1p2,
                         NULL,
        (      void ***) &tmp_ppts_field,
                         &request_elt_pts);

  PDM_part_to_part_iexch_wait(ptp, request_elt_pts);
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    free(pface_field1_p1p2[ipart]);
  }
  free(pface_field1_p1p2);



  double **ppts_field1 = malloc(sizeof(double *) * n_part_cloud);
  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    ppts_field1[ipart] = malloc(sizeof(double) * pn_pts[ipart]);

    for (int i = 0; i < pn_pts[ipart]; i++) {
      ppts_field1[ipart][i] = 0.;
    }

    for (int i = 0; i < n_located[ipart]; i++) {
      int pt_id = located[ipart][i] - 1;

      ppts_field1[ipart][pt_id] = tmp_ppts_field[ipart][i];
    }

    free(tmp_ppts_field[ipart]);
  }
  free(tmp_ppts_field);



  /*
   *  Exchange mesh <- point cloud
   */
  double **ppts_field2 = malloc(sizeof(double *) * n_part_cloud);
  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    ppts_field2[ipart] = malloc(sizeof(double) * pn_pts[ipart]);

    for (int i = 0; i < pn_pts[ipart]; i++) {
      ppts_field2[ipart][i] = cos(3*ppts_coord[ipart][3*i]);
    }
  }

  int request_pts_elt;
  double **tmp_pface_field = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(double),
                                 NULL,
                (const void  **) ppts_field2,
                                 NULL,
                (      void ***) &tmp_pface_field,
                                 &request_pts_elt);

  PDM_part_to_part_reverse_iexch_wait(ptp, request_pts_elt);

  double **pface_field2 = malloc(sizeof(double *) * n_part_mesh);
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    pface_field2[ipart] = malloc(sizeof(double) * pn_face[ipart]);

    for (int i = 0; i < pn_face[ipart]; i++) {
      pface_field2[ipart][i] = 0.;

      for (int idx = pelt_pts_idx[ipart][i]; idx < pelt_pts_idx[ipart][i+1]; idx++) {
        pface_field2[ipart][i] += tmp_pface_field[ipart][idx];
      }

      if (pelt_pts_idx[ipart][i+1] > pelt_pts_idx[ipart][i]) {
        pface_field2[ipart][i] /= (double) (pelt_pts_idx[ipart][i+1] - pelt_pts_idx[ipart][i]);
      }
    }

    free(tmp_pface_field[ipart]);
  }
  free(tmp_pface_field);



  /*
   *  Exchange mesh (node-centered) -> point cloud
   */
  double **pvtx_field3      = malloc(sizeof(double *) * n_part_mesh);
  double **pface_pts_field3 = malloc(sizeof(double *) * n_part_mesh);
  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    pvtx_field3     [ipart] = malloc(sizeof(double) * pn_vtx[ipart]);
    for (int i = 0; i < pn_vtx[ipart]; i++) {
      pvtx_field3[ipart][i] = cos(6*pvtx_coord[ipart][3*i+2]);
    }

    pface_pts_field3[ipart] = malloc(sizeof(double) * pelt_pts_idx[ipart][pn_face[ipart]]);
    for (int i = 0; i < pn_face[ipart]; i++) {

      int *fv = pface_vtx[ipart] + pface_vtx_idx[ipart][i];

      for (int idx_pt = pelt_pts_idx[ipart][i]; idx_pt < pelt_pts_idx[ipart][i+1]; idx_pt++) {
        pface_pts_field3[ipart][idx_pt] = 0.;

        int idx_vtx = 0;
        for (int idx_w = pelt_pts_weight_idx[ipart][idx_pt]; idx_w < pelt_pts_weight_idx[ipart][idx_pt+1]; idx_w++) {
          pface_pts_field3[ipart][idx_pt] += pelt_pts_weight[ipart][idx_w] * pvtx_field3[ipart][fv[idx_vtx++]-1];
        }
      }
    }
  }

  int request_elt_pts3;
  double **tmp_ppts_field3 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) pface_pts_field3,
                         NULL,
        (      void ***) &tmp_ppts_field3,
                         &request_elt_pts3);

  PDM_part_to_part_iexch_wait(ptp, request_elt_pts3);

  double **ppts_field3 = malloc(sizeof(double *) * n_part_cloud);
  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    ppts_field3[ipart] = malloc(sizeof(double) * pn_pts[ipart]);

    for (int i = 0; i < pn_pts[ipart]; i++) {
      ppts_field3[ipart][i] = 0.;
    }

    for (int i = 0; i < n_located[ipart]; i++) {
      int pt_id = located[ipart][i] - 1;

      ppts_field3[ipart][pt_id] = tmp_ppts_field3[ipart][i];
    }

    free(tmp_ppts_field3[ipart]);
  }
  free(tmp_ppts_field3);



  /* Visu */
  if (post) {
    const char    *elt_field_name[] = {"field1", "field2"};
    double **elt_field_value[2] = {pface_field1, pface_field2};
    const char    *vtx_field_name[] = {"field3"};
    double **vtx_field_value[1] = {pvtx_field3};
    _dump_mesh("proj_pt_cloud_surf__mesh",
               comm,
               n_part_mesh,
               pn_vtx,
               pvtx_coord,
               pvtx_ln_to_gn,
               pn_face,
               pface_vtx_idx,
               pface_vtx,
               pface_ln_to_gn,
               2,
               elt_field_name,
               elt_field_value,
               1,
               vtx_field_name,
               vtx_field_value);

    double **ppts_located = malloc(sizeof(double *) * n_part_cloud);
    for (int ipart = 0; ipart < n_part_cloud; ipart++) {
      ppts_located[ipart] = malloc(sizeof(double) * pn_pts[ipart]);
      for (int i = 0; i < pn_pts[ipart]; i++) {
        ppts_located[ipart][i] = 0.;
      }

      for (int i = 0; i < n_located[ipart]; i++) {
        ppts_located[ipart][located[ipart][i]-1] = 1.;
      }
    }

    const char    *field_name[] = {"field1", "field2", "field3", "located"};
    double **pts_field_value[4] = {ppts_field1, ppts_field2, ppts_field3, ppts_located};
    if (0) {
      _dump_point_cloud("proj_pt_cloud_surf__cloud",
                        comm,
                        n_part_cloud,
                        pn_pts,
                        ppts_coord,
                        ppts_ln_to_gn,
                        3,
                        field_name,
                        pts_field_value);
    }
    else {
      _dump_mesh("proj_pt_cloud_surf__cloud",
                 comm,
                 n_part_cloud,
                 pn_pts,
                 ppts_coord,
                 ppts_ln_to_gn,
                 pn_face2,
                 pface_vtx_idx2,
                 pface_vtx2,
                 pface_ln_to_gn2,
                 0,
                 NULL,
                 NULL,
                 4,
                 field_name,
                 pts_field_value);
    }

    for (int ipart = 0; ipart < n_part_cloud; ipart++) {
      free(ppts_located[ipart]);
    }
    free(ppts_located);
  }

  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    free(pface_field1    [ipart]);
    free(pface_field2    [ipart]);
    free(pvtx_field3     [ipart]);
    free(pface_pts_field3[ipart]);
  }
  free(pface_field1);
  free(pface_field2);
  free(pvtx_field3);
  free(pface_pts_field3);

  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    free(ppts_field1[ipart]);
    free(ppts_field2[ipart]);
    free(ppts_field3[ipart]);
  }
  free(ppts_field1);
  free(ppts_field2);
  free(ppts_field3);


  /* Free memory */
  PDM_part_to_part_free(ptp);
  PDM_mesh_location_free(mesh_loc);
  PDM_part_mesh_nodal_free(mesh_nodal);

  for (int ipart = 0; ipart < n_part_mesh; ipart++) {
    free(pvtx_coord    [ipart]);
    free(pvtx_ln_to_gn [ipart]);
    free(pface_vtx_idx [ipart]);
    free(pface_vtx     [ipart]);
    free(pface_ln_to_gn[ipart]);
  }
  free(pvtx_coord);
  free(pvtx_ln_to_gn);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pface_ln_to_gn);
  free(pn_vtx);
  free(pn_face);

  free(pelt_pts_idx);
  free(pelt_pts_gnum);
  free(pelt_pts_coord);
  free(pelt_pts_uvw);
  free(pelt_pts_weight_idx);
  free(pelt_pts_weight);
  free(pelt_pts_dist2);
  free(pelt_pts_proj_coord);

  for (int ipart = 0; ipart < n_part_cloud; ipart++) {
    free(ppts_coord    [ipart]);
    free(ppts_ln_to_gn [ipart]);
    free(pface_vtx_idx2 [ipart]);
    free(pface_vtx2     [ipart]);
    free(pface_ln_to_gn2[ipart]);
  }
  free(ppts_coord);
  free(ppts_ln_to_gn);
  free(pn_pts);
  free(pface_vtx_idx2);
  free(pface_vtx2);
  free(pface_ln_to_gn2);
  free(pn_face2);


  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
