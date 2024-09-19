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
#include "pdm_sort.h"
#include "pdm_dcube_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_part_to_part.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"


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
           double        *separation,
           double        *tolerance,
           double        *marge,
           int           *n_part,
           PDM_g_num_t   *n_pts,
           int           *post,
           int           *part_method,
           PDM_mesh_location_method_t *loc_method)
{
  int i = 1;

  PDM_UNUSED (post);

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
    else if (strcmp(argv[i], "-sep") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation = atof(argv[i]);
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
        *n_pts = atoi(argv[i]);
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
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static PDM_part_t *
_cube_mesh
(
 const int              n_part,
 const PDM_part_split_t part_method,
 const PDM_g_num_t      n_vtx_seg,
 const double           xmin,
 const double           ymin,
 const double           zmin,
 const double           length
 )
{
  PDM_dcube_t *dcube = PDM_dcube_gen_init(PDM_MPI_COMM_WORLD,
                                          n_vtx_seg,
                                          length,
                                          xmin,
                                          ymin,
                                          zmin,
                                          PDM_OWNERSHIP_KEEP);

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


  PDM_dcube_gen_dim_get (dcube,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtx_l,
                         &dface_group_l);

  PDM_dcube_gen_data_get (dcube,
                          &dface_cell,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &dvtx_coord,
                          &dface_group_idx,
                          &dface_group);

  /*
   *  Create mesh partitiions
   */
  // int ppart_id = 0;
  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_t *ppart = PDM_part_create (PDM_MPI_COMM_WORLD,
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

  PDM_dcube_gen_free(dcube);

  return ppart;
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

  PDM_g_num_t n_vtx_seg = 3;
  double      length    = 1.;
  double      separation = 0.5;
  double      tolerance = 1e-3;
  double      marge     = 0.;
  int         n_part    = 1;
  int         post      = 0;
  PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;

  PDM_g_num_t n_pts = 10;
  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &length,
              &separation,
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
  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);


  const double xmin = 0;
  const double ymin = 0;
  const double zmin = 0;

  //double step = length / (n_vtx_seg - 1);
  double delta = separation * length;

  const double xmin2 = length - delta ;
  const double ymin2 = 0.;
  const double zmin2 = 0.;
  double length2 = length;

  if (i_rank == 0) {
  printf("xmin1, ymin1, zmin1 : %12.5e %12.5e %12.5e\n", xmin, ymin,zmin);
  printf("xmax1, ymax1, zmax1 : %12.5e %12.5e %12.5e\n", xmin+length, ymin+length,zmin+length);
  PDM_g_num_t n_elt = n_vtx_seg * n_vtx_seg * n_vtx_seg;
  printf("n_elt : "PDM_FMT_G_NUM"\n", n_elt);

  printf("xmin2, ymin2, zmin2 : %12.5e %12.5e %12.5e\n", xmin2, ymin2,zmin2);
  printf("xmax2, ymax2, zmax2 : %12.5e %12.5e %12.5e\n", xmin2+length2, ymin2+length2,zmin2+length2);
  }
  /*
   *  Source cube
   */
  PDM_part_t *ppart_src = _cube_mesh (n_part,
                                      part_method,
                                      n_vtx_seg,
                                      xmin,
                                      ymin,
                                      zmin,
                                      length);

  /*
   *  Target cube
   */
  PDM_part_t *ppart_tgt = _cube_mesh (n_part,
                                      part_method,
                                      2*n_vtx_seg,
 //                             n_vtx_seg,
                                      xmin2,
                                      ymin2,
                                      zmin2,
                                      length2);

  /*
   *  Mesh location structure initialization
   */
  PDM_mesh_location_t *id_loc1 = PDM_mesh_location_create(1,
                                                          PDM_MPI_COMM_WORLD,
                                                          PDM_OWNERSHIP_KEEP);

  PDM_mesh_location_t *id_loc2 = PDM_mesh_location_create(1,
                                                          PDM_MPI_COMM_WORLD,
                                                          PDM_OWNERSHIP_KEEP);

  /* Set target point cloud */
  PDM_mesh_location_n_part_cloud_set (id_loc1,
                                      0,
                                      n_part);

  PDM_mesh_location_n_part_cloud_set (id_loc2,
                                      0,
                                      n_part);
  double **cell_volume1 = malloc (sizeof(double *) * n_part);
  double **cell_center1 = malloc (sizeof(double *) * n_part);

  double **cell_volume2 = malloc (sizeof(double *) * n_part);
  double **cell_center2 = malloc (sizeof(double *) * n_part);

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

    PDM_part_part_dim_get (ppart_tgt,
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

    PDM_part_part_val_get (ppart_tgt,
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

    const int is_oriented = 0;
    cell_volume2[ipart] = malloc(sizeof(double) * n_cell);
    cell_center2[ipart] = malloc(sizeof(double) * 3 * n_cell);

    PDM_geom_elem_polyhedra_properties (is_oriented,
                                        n_cell,
                                        n_face,
                                        face_vtx_idx,
                                        face_vtx,
                                        cell_face_idx,
                                        cell_face,
                                        n_vtx,
                                        vtx,
                                        cell_volume2[ipart],
                                        cell_center2[ipart],
                                        NULL,
                                        NULL);

    PDM_mesh_location_cloud_set (id_loc1,
                                 0,
                                 ipart,
                                 n_cell,
                                 cell_center2[ipart],
                                 cell_ln_to_gn);
  }



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

    PDM_part_part_dim_get (ppart_src,
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

    PDM_part_part_val_get (ppart_src,
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

    const int is_oriented = 0;
    cell_volume1[ipart] = malloc(sizeof(double) * n_cell);
    cell_center1[ipart] = malloc(sizeof(double) * 3 * n_cell);

    PDM_geom_elem_polyhedra_properties (is_oriented,
                                        n_cell,
                                        n_face,
                                        face_vtx_idx,
                                        face_vtx,
                                        cell_face_idx,
                                        cell_face,
                                        n_vtx,
                                        vtx,
                                        cell_volume1[ipart],
                                        cell_center1[ipart],
                                        NULL,
                                        NULL);

    PDM_mesh_location_cloud_set (id_loc2,
                                 0,
                                 ipart,
                                 n_cell,
                                 cell_center1[ipart],
                                 cell_ln_to_gn);
  }




  /* Set source mesh */
  PDM_mesh_location_mesh_n_part_set (id_loc1,
                                          n_part);

  PDM_mesh_location_mesh_n_part_set (id_loc2,
                                          n_part);

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

    PDM_part_part_dim_get (ppart_src,
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

    PDM_part_part_val_get (ppart_src,
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

    PDM_mesh_location_part_set (id_loc1,
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

    PDM_part_part_dim_get (ppart_tgt,
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

    PDM_part_part_val_get (ppart_tgt,
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

    PDM_mesh_location_part_set (id_loc2,
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
  PDM_mesh_location_tolerance_set (id_loc1,
                                   tolerance);

  PDM_mesh_location_method_set (id_loc1,
                                loc_method);

  PDM_mesh_location_tolerance_set (id_loc2,
                                   tolerance);

  PDM_mesh_location_method_set (id_loc2,
                                loc_method);


  /*
   *  Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  PDM_mesh_location_compute (id_loc1);

  PDM_mesh_location_dump_times (id_loc1);

  PDM_mesh_location_compute (id_loc2);

  PDM_mesh_location_dump_times (id_loc2);

  int         **elt_pts_inside_idx      = malloc (sizeof(int         *) * n_part);
  PDM_g_num_t **location                = malloc (sizeof(PDM_g_num_t *) * n_part);
  int         **location_idx            = malloc (sizeof(int         *) * n_part);
  PDM_g_num_t **points_gnum             = malloc (sizeof(PDM_g_num_t *) * n_part);
  double      **points_coords           = malloc (sizeof(double      *) * n_part);
  double      **points_uvw              = malloc (sizeof(double      *) * n_part);
  int         **points_weights_idx      = malloc (sizeof(int         *) * n_part);
  double      **points_weights          = malloc (sizeof(double      *) * n_part);
  double      **points_dist2            = malloc (sizeof(double      *) * n_part);
  double      **points_projected_coords = malloc (sizeof(double      *) * n_part);
  PDM_g_num_t **gnum_elt1               = malloc (sizeof(PDM_g_num_t *) * n_part);
  int          *n_elt1                  = malloc (sizeof(int          ) * n_part);
  PDM_g_num_t **gnum_elt2               = malloc (sizeof(PDM_g_num_t *) * n_part);
  int          *n_elt2                  = malloc (sizeof(int          ) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {

    int n_face_tgt;
    int n_face_part_bound_tgt;
    int n_vtx_tgt;
    int n_proc_tgt;
    int n_t_part_tgt;
    int s_cell_face_tgt;
    int s_face_vtx_tgt;
    int s_face_group_tgt;
    int n_edge_group2_tgt;

    PDM_part_part_dim_get (ppart_tgt,
                           ipart,
                           &(n_elt2[ipart]),
                           &n_face_tgt,
                           &n_face_part_bound_tgt,
                           &n_vtx_tgt,
                           &n_proc_tgt,
                           &n_t_part_tgt,
                           &s_cell_face_tgt,
                           &s_face_vtx_tgt,
                           &s_face_group_tgt,
                           &n_edge_group2_tgt);

    location_idx[ipart] = malloc (sizeof(int) * (n_elt2[ipart]+1));
    for (int i = 0; i < n_elt2[ipart]+1; i++) {
      location_idx[ipart][i] = 0;
    }
    int n_located = PDM_mesh_location_n_located_get (id_loc1,
                                                     0,//i_point_cloud,
                                                     ipart);

    int *located = PDM_mesh_location_located_get (id_loc1,
                                                  0,//i_point_cloud,
                                                  ipart);
    for (int i = 0; i < n_located; i++) {
      location_idx[ipart][located[i]-1+1] = 1;
    }
    for (int i = 0; i < n_elt2[ipart]; i++) {
      location_idx[ipart][i+1] += location_idx[ipart][i]; 
    }

    int n_unlocated = PDM_mesh_location_n_unlocated_get (id_loc1,
                                                         0,//i_point_cloud,
                                                         ipart);

    int *unlocated = PDM_mesh_location_unlocated_get (id_loc1,
                                                      0,//i_point_cloud,
                                                      ipart);

    if (1 == 0) {
      printf ("located :");
      for (int i = 0; i < n_located; i++) {
        printf (" %d", located[i]);
      }
      printf("\n");

      printf ("unlocated :");
      for (int i = 0; i < n_unlocated; i++) {
        printf (" %d", unlocated[i]);
      }
      printf("\n");
    }

    int         *cell_tag_tgt;
    int         *cell_face_idx_tgt;
    int         *cell_face_tgt;
    int         *face_tag_tgt;
    int         *face_cell_tgt;
    int         *face_vtx_idx_tgt;
    int         *face_vtx_tgt;
    PDM_g_num_t *face_ln_to_gn_tgt;
    int         *face_part_boundProcIdx_tgt;
    int         *face_part_boundPartIdx_tgt;
    int         *face_part_bound_tgt;
    int         *vtx_tag_tgt;
    double      *vtx_tgt;
    PDM_g_num_t *vtx_ln_to_gn_tgt;
    int         *face_group_idx_tgt;
    int         *face_group_tgt;
    PDM_g_num_t *face_group_ln_to_gn_tgt;

    PDM_part_part_val_get (ppart_tgt,
                           ipart,
                           &cell_tag_tgt,
                           &cell_face_idx_tgt,
                           &cell_face_tgt,
                           &(gnum_elt2[ipart]),
                           &face_tag_tgt,
                           &face_cell_tgt,
                           &face_vtx_idx_tgt,
                           &face_vtx_tgt,
                           &face_ln_to_gn_tgt,
                           &face_part_boundProcIdx_tgt,
                           &face_part_boundPartIdx_tgt,
                           &face_part_bound_tgt,
                           &vtx_tag_tgt,
                           &vtx_tgt,
                           &vtx_ln_to_gn_tgt,
                           &face_group_idx_tgt,
                           &face_group_tgt,
                           &face_group_ln_to_gn_tgt);



    int n_face_src;
    int n_face_part_bound_src;
    int n_vtx_src;
    int n_proc_src;
    int n_t_part_src;
    int s_cell_face_src;
    int s_face_vtx_src;
    int s_face_group_src;
    int n_edge_group2_src;

    PDM_part_part_dim_get (ppart_src,
                           ipart,
                           &(n_elt1[ipart]),
                           &n_face_src,
                           &n_face_part_bound_src,
                           &n_vtx_src,
                           &n_proc_src,
                           &n_t_part_src,
                           &s_cell_face_src,
                           &s_face_vtx_src,
                           &s_face_group_src,
                           &n_edge_group2_src);

    int         *cell_tag_src;
    int         *cell_face_idx_src;
    int         *cell_face_src;
    int         *face_tag_src;
    int         *face_cell_src;
    int         *face_vtx_idx_src;
    int         *face_vtx_src;
    PDM_g_num_t *face_ln_to_gn_src;
    int         *face_part_boundProcIdx_src;
    int         *face_part_boundPartIdx_src;
    int         *face_part_bound_src;
    int         *vtx_tag_src;
    double      *vtx_src;
    PDM_g_num_t *vtx_ln_to_gn_src;
    int         *face_group_idx_src;
    int         *face_group_src;
    PDM_g_num_t *face_group_ln_to_gn_src;

    PDM_part_part_val_get (ppart_src,
                           ipart,
                           &cell_tag_src,
                           &cell_face_idx_src,
                           &cell_face_src,
                           &(gnum_elt1[ipart]),
                           &face_tag_src,
                           &face_cell_src,
                           &face_vtx_idx_src,
                           &face_vtx_src,
                           &face_ln_to_gn_src,
                           &face_part_boundProcIdx_src,
                           &face_part_boundPartIdx_src,
                           &face_part_bound_src,
                           &vtx_tag_src,
                           &vtx_src,
                           &vtx_ln_to_gn_src,
                           &face_group_idx_src,
                           &face_group_src,
                           &face_group_ln_to_gn_src);



    PDM_mesh_location_points_in_elt_get (id_loc1,
                                        ipart,
                                        0,
                                       &(elt_pts_inside_idx[ipart]),
                                       &(points_gnum[ipart]),
                                       &(points_coords[ipart]),
                                       &(points_uvw[ipart]),
                                       &(points_weights_idx[ipart]),
                                       &(points_weights[ipart]),
                                       &(points_dist2[ipart]),
                                       &(points_projected_coords[ipart]));


    if (1 == 0) {
      printf ("elt_pts_inside :\n");
      for (int i = 0; i < n_elt1[ipart]; i++) {
        printf ("%d "PDM_FMT_G_NUM": ", i, gnum_elt1[ipart][i]);
        for (int j = elt_pts_inside_idx[ipart][i]; j < elt_pts_inside_idx[ipart][i+1]; j++) {
          printf (" "PDM_FMT_G_NUM"", points_gnum[ipart][j]);
        }
        printf ("\n");
      }
      printf("\n");
    }

    double              *dist2;
    double              *projected_coord;

    PDM_mesh_location_point_location_get (id_loc1,
                                          0,
                                          ipart,
                                          &(location[ipart]),
                                          &dist2,
                                          &projected_coord);

    //assert (n_located == 0 && n_unlocated == n_cell);
    free (cell_center2[ipart]);
    free (cell_volume2[ipart]);
    free (cell_center1[ipart]);
    free (cell_volume1[ipart]);

  }


  PDM_part_to_part_t *ptp = PDM_part_to_part_create ((const PDM_g_num_t**) gnum_elt1,
                                                                         n_elt1,
                                                                         n_part,
                                                                         (const PDM_g_num_t**) gnum_elt2,
                                                                         n_elt2,
                                                                         n_part,
                                                                         (const int **) elt_pts_inside_idx,
                                                                         (const PDM_g_num_t **) points_gnum,
                                                                         PDM_MPI_COMM_WORLD);


  PDM_part_to_part_t *ptp2 = PDM_part_to_part_create ((const PDM_g_num_t**) gnum_elt2,
                                                                         n_elt2,
                                                                         n_part,
                                                                         (const PDM_g_num_t**) gnum_elt1,
                                                                         n_elt1,
                                                                         n_part,
                                                                         (const int **) location_idx,
                                                                         (const PDM_g_num_t **) location,
                                                                         PDM_MPI_COMM_WORLD);

  int  *n_ref_lnum2;
  int **ref_lnum2;
  PDM_part_to_part_ref_lnum2_get (ptp,
                                  &n_ref_lnum2,
                                  &ref_lnum2);


  int  *n_unref_lnum2;
  int **unref_lnum2;
  PDM_part_to_part_unref_lnum2_get (ptp,
                                    &n_unref_lnum2,
                                    &unref_lnum2);


  int         **gnum1_come_from_idx;
  PDM_g_num_t **gnum1_come_from;
  PDM_part_to_part_gnum1_come_from_get (ptp,
                                        &gnum1_come_from_idx,
                                        &gnum1_come_from);


  int  *ptp2_n_ref_lnum2;
  int **ptp2_ref_lnum2;
  PDM_part_to_part_ref_lnum2_get (ptp2,
                                  &ptp2_n_ref_lnum2,
                                  &ptp2_ref_lnum2);

  int  *ptp2_n_unref_lnum2;
  int **ptp2_unref_lnum2;
  PDM_part_to_part_unref_lnum2_get (ptp2,
                                    &ptp2_n_unref_lnum2,
                                    &ptp2_unref_lnum2);


  int         **ptp2_gnum1_come_from_idx;
  PDM_g_num_t **ptp2_gnum1_come_from;
  PDM_part_to_part_gnum1_come_from_get (ptp2,
                                        &ptp2_gnum1_come_from_idx,
                                        &ptp2_gnum1_come_from);


  int send_request = -1;

  PDM_g_num_t **gnum1_gnum2_data = malloc (sizeof(PDM_g_num_t *) * n_part);
  for (int i = 0; i < n_part; i++) {
    gnum1_gnum2_data[i] = malloc (sizeof(PDM_g_num_t) * elt_pts_inside_idx[i][n_elt1[i]]);
    for (int j = 0; j < n_elt1[i]; j++) {
      for (int k = elt_pts_inside_idx[i][j]; k < elt_pts_inside_idx[i][j+1]; k++) {
        gnum1_gnum2_data[i][k] = gnum_elt1[i][j];
      }
    }
  }

  PDM_part_to_part_issend (ptp,
                           sizeof (PDM_g_num_t),
                           1,
         (const void **)  gnum1_gnum2_data,
                           100,
                           &send_request);



  int recv_request = -1;
  PDM_g_num_t **gnum_elt1_recv = malloc (sizeof(PDM_g_num_t*)  * n_part);
  for (int i = 0; i < n_part; i++) {
    gnum_elt1_recv[i] = malloc (sizeof(PDM_g_num_t)  * gnum1_come_from_idx[i][n_ref_lnum2[i]]);
  }

  PDM_part_to_part_irecv (ptp,
                          sizeof (PDM_g_num_t),
                          1,
                (void **) gnum_elt1_recv,
                          100,
                          &recv_request);



  PDM_part_to_part_issend_wait (ptp, send_request);

  PDM_part_to_part_irecv_wait (ptp, recv_request);

  for (int i = 0; i < n_part; i++) {
    free (gnum1_gnum2_data[i]);
  }
  free (gnum1_gnum2_data);


  PDM_g_num_t **ptp2_s_data = malloc (sizeof(PDM_g_num_t *) * n_part);
  for (int i = 0; i < n_part; i++) {
    ptp2_s_data[i] = malloc (sizeof(PDM_g_num_t) * location_idx[i][n_elt2[i]]);
    for (int j = 0; j < n_elt2[i]; j++) {
      for (int k = location_idx[i][j]; k < location_idx[i][j+1]; k++) {
        ptp2_s_data[i][k] = gnum_elt2[i][j];
      }
    }
  }

  send_request = -1;
  PDM_part_to_part_issend (ptp2,
                           sizeof (PDM_g_num_t),
                           1,
           (const void **) ptp2_s_data,
                           100,
                           &send_request);

  recv_request = -1;
  PDM_g_num_t **gnum_elt2_recv = malloc (sizeof(PDM_g_num_t*)  * n_part);
  for (int i = 0; i < n_part; i++) {
    gnum_elt2_recv[i] = malloc (sizeof(PDM_g_num_t)  * ptp2_gnum1_come_from_idx[i][ptp2_n_ref_lnum2[i]]);
  }

  PDM_part_to_part_irecv (ptp2,
                          sizeof (PDM_g_num_t),
                          1,
                (void **) gnum_elt2_recv,
                          100,
                          &recv_request);



  PDM_part_to_part_issend_wait (ptp2, send_request);

  PDM_part_to_part_irecv_wait (ptp2, recv_request);

  for (int i = 0; i < n_part; i++) {
    free (ptp2_s_data[i]);
  }
  free (ptp2_s_data);


  for (int i = 0; i < n_part; i++) {
    int n_located = PDM_mesh_location_n_located_get (id_loc1,
                                                     0,//i_point_cloud,
                                                     i);

    int *located = PDM_mesh_location_located_get (id_loc1,
                                                  0,//i_point_cloud,
                                                  i);
    assert (n_located == n_ref_lnum2[i]);

    if (1 == 0) {
      printf ("located :");
      for (int j = 0; j < n_located; j++) {
        printf(" %d", located[j]);
      }
      printf ("\n");

      printf ("ref_lnum2 :");
      for (int j = 0; j < n_ref_lnum2[i]; j++) {
        for (int k = gnum1_come_from_idx[i][j] ; k < gnum1_come_from_idx[i][j+1]; k++) {
          printf(" "PDM_FMT_G_NUM"", gnum1_come_from[i][k]);
        }
        printf ("\n");
      }
    }

 /*   PDM_g_num_t         *location;
    double              *dist2;
    double              *projected_coord;

    PDM_mesh_location_point_location_get (id_loc1,
                                          0,
                                          i,
                                          &location,
                                          &dist2,
                                          &projected_coord);
*/
    int n_face_tgt;
    int n_face_part_bound_tgt;
    int n_vtx_tgt;
    int n_proc_tgt;
    int n_t_part_tgt;
    int s_cell_face_tgt;
    int s_face_vtx_tgt;
    int s_face_group_tgt;
    int n_edge_group2_tgt;

    PDM_part_part_dim_get (ppart_tgt,
                           i,
                           &(n_elt2[i]),
                           &n_face_tgt,
                           &n_face_part_bound_tgt,
                           &n_vtx_tgt,
                           &n_proc_tgt,
                           &n_t_part_tgt,
                           &s_cell_face_tgt,
                           &s_face_vtx_tgt,
                           &s_face_group_tgt,
                           &n_edge_group2_tgt);

    int         *cell_tag_tgt;
    int         *cell_face_idx_tgt;
    int         *cell_face_tgt;
    int         *face_tag_tgt;
    int         *face_cell_tgt;
    int         *face_vtx_idx_tgt;
    int         *face_vtx_tgt;
    PDM_g_num_t *face_ln_to_gn_tgt;
    int         *face_part_boundProcIdx_tgt;
    int         *face_part_boundPartIdx_tgt;
    int         *face_part_bound_tgt;
    int         *vtx_tag_tgt;
    double      *vtx_tgt;
    PDM_g_num_t *vtx_ln_to_gn_tgt;
    int         *face_group_idx_tgt;
    int         *face_group_tgt;
    PDM_g_num_t *face_group_ln_to_gn_tgt;

    PDM_part_part_val_get (ppart_tgt,
                           i,
                           &cell_tag_tgt,
                           &cell_face_idx_tgt,
                           &cell_face_tgt,
                           &(gnum_elt2[i]),
                           &face_tag_tgt,
                           &face_cell_tgt,
                           &face_vtx_idx_tgt,
                           &face_vtx_tgt,
                           &face_ln_to_gn_tgt,
                           &face_part_boundProcIdx_tgt,
                           &face_part_boundPartIdx_tgt,
                           &face_part_bound_tgt,
                           &vtx_tag_tgt,
                           &vtx_tgt,
                           &vtx_ln_to_gn_tgt,
                           &face_group_idx_tgt,
                           &face_group_tgt,
                           &face_group_ln_to_gn_tgt);

    for (int j = 0; j < n_elt1[i]; j++) {
      int n_pts2 = elt_pts_inside_idx[i][j+1] - elt_pts_inside_idx[i][j]; 
      PDM_sort_long (points_gnum[i] + elt_pts_inside_idx[i][j], NULL, n_pts2);
    }

    if (0 == 1) {
 
      printf ("location from location :\n");
      for (int j = 0; j < n_located; j++) {
        printf(""PDM_FMT_G_NUM" : "PDM_FMT_G_NUM"", gnum_elt2[i][located[j]-1], location[i][j]);
      printf ("\n");
      }
 
      printf ("location from exchange  :\n");
      for (int j = 0; j < n_ref_lnum2[i]; j++) {
        printf(""PDM_FMT_G_NUM" :", gnum_elt2[i][ref_lnum2[i][j]-1]);
        for (int k = gnum1_come_from_idx[i][j] ; k < gnum1_come_from_idx[i][j+1]; k++) {
          printf(" "PDM_FMT_G_NUM"", gnum_elt1_recv[i][k]);
        }  
        printf ("\n");
      }

      printf ("elt_pts_inside from location :\n");
      for (int j = 0; j < n_elt1[i]; j++) {
        printf ("%d "PDM_FMT_G_NUM": ", j, gnum_elt1[i][j]);
        int n_pts2 = elt_pts_inside_idx[i][j+1] - elt_pts_inside_idx[i][j]; 
        PDM_sort_long (points_gnum[i] + elt_pts_inside_idx[i][j], NULL, n_pts2);

        for (int k = elt_pts_inside_idx[i][j]; k < elt_pts_inside_idx[i][j+1]; k++) {
          printf (" "PDM_FMT_G_NUM"", points_gnum[i][k]);
        }
        printf ("\n");
      }
      printf("\n");

      printf ("elt_pts_inside from exchange :\n");
      for (int j = 0; j < ptp2_n_ref_lnum2[i]; j++) {
        printf(""PDM_FMT_G_NUM" :", gnum_elt1[i][ptp2_ref_lnum2[i][j]-1]);
        for (int k = ptp2_gnum1_come_from_idx[i][j] ; k < ptp2_gnum1_come_from_idx[i][j+1]; k++) {
          printf(" "PDM_FMT_G_NUM"", gnum_elt2_recv[i][k]);
        }  
        printf ("\n");
      }
    }

    for (int j = 0; j < n_ref_lnum2[i]; j++) {
      for (int k = gnum1_come_from_idx[i][j] ; k < gnum1_come_from_idx[i][j+1]; k++) {
        assert(gnum_elt1_recv[i][k] == location[i][j]);
      }  
    }

    for (int j = 0; j < ptp2_n_ref_lnum2[i]; j++) {
      int ielt = ptp2_ref_lnum2[i][j]-1;
      int n1 = ptp2_gnum1_come_from_idx[i][j+1] - ptp2_gnum1_come_from_idx[i][j];
      int n2 = elt_pts_inside_idx[i][ielt+1] - elt_pts_inside_idx[i][ielt];  
      assert(n1 == n2);
      for (int k = ptp2_gnum1_come_from_idx[i][j], k1 =  elt_pts_inside_idx[i][ielt] ; k < ptp2_gnum1_come_from_idx[i][j+1]; k++, k1++) {
        PDM_g_num_t val1 = gnum_elt2_recv[i][k];
        PDM_g_num_t val2 = points_gnum[i][k1];
        assert( val1 == val2); 
      }  
    }
  }

  for (int i = 0; i < n_part; i++) {
    free (gnum_elt1_recv[i]);
    free (gnum_elt2_recv[i]);
    free (location_idx[i]);
  }

  free (gnum_elt1_recv);
  free (gnum_elt2_recv);

  free (cell_center1);
  free (cell_volume1);
  free (cell_center2);
  free (cell_volume2);
  free (location_idx);
  free (location);

  free (elt_pts_inside_idx);
  free (points_gnum);
  free (points_coords); 
  free (points_uvw); 
  free (points_weights_idx);
  free (points_weights);
  free (points_dist2); 
  free (points_projected_coords);

  free (gnum_elt1);
  free (n_elt1);
  free (gnum_elt2);
  free (n_elt2);

  PDM_mesh_location_free (id_loc1);

  PDM_mesh_location_free (id_loc2);

  PDM_part_to_part_free (ptp);
  PDM_part_to_part_free (ptp2);

  PDM_part_free (ppart_src);
  PDM_part_free (ppart_tgt);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
