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
#include "pdm_multipart.h"
#include "pdm_dcube_gen.h"
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_part_to_block.h"

#include "pdm_part_extension.h"
#include "pdm_vtk.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_point_cloud_gen.h"

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
     "  -pt-scotch       Call PT-Scotch.\n\n"
     "  -hilbert         Call Hilbert.\n\n"
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
_read_args(int                          argc,
           char                       **argv,
           PDM_g_num_t                 *n_vtx_seg,
           double                      *length,
           double                      *separation_x,
           double                      *separation_y,
           double                      *separation_z,
           int                         *deform,
           double                      *tolerance,
           double                      *marge,
           int                         *n_part,
           PDM_g_num_t                 *n_pts,
           int                         *post,
           int                         *part_method,
           PDM_mesh_location_method_t  *loc_method,
           int                         *order,
           PDM_Mesh_nodal_elt_t        *elt_type)
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
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepy") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_y = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_z = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-def") == 0) {
      *deform = 1;
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
    else if (strcmp(argv[i], "-doctree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DOCTREE;
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-order") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *order = atoi(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

// static double R[3][3] =
// {
//   {-0.14547275709949994,  0.8415293589391187 , -0.5202557207618055 },
//   { 0.9893622576902102 ,  0.12373586628506748, -0.07649678720582984},
//   { 0.                 , -0.5258495730132333 , -0.8505775840931856 }
// };

// static void _rotate
// (
//  const int     n_pts,
//        double *coord
//  )
// {
//   for (int i = 0; i < n_pts; i++) {
//     double x = coord[3*i];
//     double y = coord[3*i+1];
//     double z = coord[3*i+2];

//     for (int j = 0; j < 3; j++) {
//       coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
//     }
//   }
// }

// static void
// _unrotate
// (
//  const int     n_pts,
//        double *coord
//  )
// {
//   for (int i = 0 ; i < n_pts ; i++) {
//     double x = coord[3 * i];
//     double y = coord[3 * i + 1];
//     double z = coord[3 * i + 2];

//     for (int j = 0 ; j < 3 ; j++) {
//       coord[3 * i + j] = R[0][j] * x + R[1][j] * y + R[2][j] * z;
//     }
//   }
// }

static void
_eval_deformation
(
 const int     order,
 const double  x,
 const double  y,
 const double  z,
       double *dx,
       double *dy,
       double *dz
 )
{
  *dx = 0.4*y;
  *dy = 0.4*z;
  *dz = 0.4*x;

  if (order == 1) {
    return;
  }

  *dx += 0.2*y*y;
  *dy += 0.2*z*z;
  *dz += 0.2*x*x;

  if (order == 2) {
    return;
  }

  *dx += 0.1*y*y*y;
  *dy += 0.1*z*z*z;
  *dz += 0.1*x*x*x;

  return;
}


static PDM_part_mesh_nodal_t *
_gen_mesh
(
 const PDM_MPI_Comm            comm,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 const PDM_g_num_t             n_vtx_seg,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  length,
 const int                     deform,
 const PDM_Mesh_nodal_elt_t    elt_type,
 const int                     order,
       PDM_multipart_t       **mpart
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_domain = 1;
  *mpart = PDM_multipart_create(n_domain,
                                &n_part,
                                PDM_FALSE,
                                part_method,
                                PDM_PART_SIZE_HOMOGENEOUS,
                                NULL,
                                comm,
                                PDM_OWNERSHIP_KEEP);


  PDM_multipart_set_reordering_options(*mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  assert(elt_type != PDM_MESH_NODAL_POLY_3D);

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         xmin,
                                                         ymin,
                                                         zmin,
                                                         elt_type,
                                                         order,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_build(dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  if (deform) {
    int order_deformation = PDM_MAX(1, (int) ceil (sqrt(order)));
    const PDM_g_num_t *distrib_vtx = PDM_DMesh_nodal_distrib_vtx_get(dmn);
    int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];
    double *dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn);
    // _rotate(dn_vtx,
    //         dvtx_coord);
    for (int i = 0; i < dn_vtx; i++) {
      double dxyz[3];
      _eval_deformation(order_deformation,
                        dvtx_coord[3*i+0],
                        dvtx_coord[3*i+1],
                        dvtx_coord[3*i+2],
                        &dxyz[0],
                        &dxyz[1],
                        &dxyz[2]);
      for (int j = 0; j < 3; j++) {
        dvtx_coord[3*i+j] += dxyz[j];
      }
    }
  }

  /*
   * Split mesh
   */
  PDM_multipart_dmesh_nodal_set(*mpart, 0, dmn);
  PDM_multipart_compute(*mpart);

  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(*mpart, 0, &pmn, PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_free(dcube);

  return pmn;
}


// static inline double
// _eval_field
// (
//  double *xyz
//  )
// {
//   return 1 + 2*xyz[0] + 3*xyz[1] + 4*xyz[2];
// }

static double
_eval_field
(
 const double x,
 const double y,
 const double z,
 const int    order
 )
{
  double f = x - y + z - 1;

  if (order == 1) {
    return f;
  }

  f += x*x - y*y + z*z - x*y + y*z - z*x;

  if (order == 2) {
    return f;
  }

  f +=
  x*x*x - y*y*y + z*z*z -
  x*x*y + x*x*z -
  y*y*x + y*y*z -
  z*z*x + z*z*y;

  return f;
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

  PDM_g_num_t n_vtx_seg           = 10;
  double      length              = 1.;
  double      separation_x        = 2.;
  double      separation_y        = 0.;
  double      separation_z        = 0.;
  int         deform              = 0;
  double      tolerance           = 1e-6;
  double      marge               = 0.;
  int         n_part              = 1;
  int         post                = 0;
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;

  PDM_g_num_t gn_pts = 10;
  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;

  int                  order    = 1;
  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_HEXA8;

  //  2 -> tria
  //  3 -> quad
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  // 11 -> tria_ho
  // 12 -> quad_ho
  // 13 -> tetra_ho
  // 14 -> pyramid_ho
  // 15 -> prism_ho
  // 16 -> hexa_ho


  int order_deformation = PDM_MAX(1, (int) ceil (sqrt(order)));
  int order_field       = PDM_MAX(1, (int) floor(sqrt(order)));

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &separation_x,
             &separation_y,
             &separation_z,
             &deform,
             &tolerance,
             &marge,
             &n_part,
             &gn_pts,
             &post,
     (int *) &part_method,
             &loc_method,
             &order,
             &elt_type);

  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 3);

  if (order > 3) {
    int *ijk = NULL;

    for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_BARHO;
         type <= PDM_MESH_NODAL_HEXAHO;
         type++) {

      if (type == PDM_MESH_NODAL_PYRAMIDHO) continue;

      ijk = PDM_vtk_lagrange_to_ijk(type, order);
      PDM_ho_ordering_user_to_ijk_add ("PDM_HO_ORDERING_VTK",
                                       type,
                                       order,
                                       PDM_Mesh_nodal_n_vtx_elt_get(type, order),
                                       ijk);
      free (ijk);
    }
  }

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  /*
   *  Source mesh
   */
  const double xmin = -0.5*length;
  const double ymin = -0.5*length;
  const double zmin = -0.5*length;

  PDM_multipart_t *src_mpart = NULL;
  PDM_part_mesh_nodal_t *src_pmn = _gen_mesh(comm,
                                             n_part,
                                             part_method,
                                             n_vtx_seg,
                                             xmin,
                                             ymin,
                                             zmin,
                                             length,
                                             deform,
                                             elt_type,
                                             order,
                                             &src_mpart);

  // if (post) {
  //   PDM_part_mesh_nodal_dump_vtk(src_pmn,
  //                                PDM_GEOMETRY_KIND_VOLUMIC,
  //                                "mesh_location_ho_src");
  // }

  /*
   *  Target cloud
   */
  int          n_pts;
  double      *pts_coord;
  PDM_g_num_t *pts_ln_to_gn;
  double h = 0.5*length/(double) (n_vtx_seg - 1);
  PDM_point_cloud_gen_random(comm,
                             0,
                             0,
                             gn_pts,
                             xmin + h,
                             ymin + h,
                             zmin + h,
                             xmin + length - h,
                             ymin + length - h,
                             zmin + length - h,
                             &n_pts,
                             &pts_coord,
                             &pts_ln_to_gn);
  // PDM_point_cloud_gen_cartesian(comm,
  //                               (int) n_vtx_seg,
  //                               (int) n_vtx_seg,
  //                               (int) n_vtx_seg,
  //                               xmin + h,
  //                               ymin + h,
  //                               zmin + h,
  //                               xmin + length - h,
  //                               ymin + length - h,
  //                               zmin + length - h,
  //                               &n_pts,
  //                               &pts_coord,
  //                               &pts_ln_to_gn);

  if (deform) {
    // _rotate(n_pts,
    //         pts_coord);
    for (int i = 0; i < n_pts; i++) {
      double dxyz[3];
      _eval_deformation(order_deformation,
                        pts_coord[3*i+0],
                        pts_coord[3*i+1],
                        pts_coord[3*i+2],
                        &dxyz[0],
                        &dxyz[1],
                        &dxyz[2]);
      for (int j = 0; j < 3; j++) {
        pts_coord[3*i+j] += dxyz[j];
      }
    }
  }


  /*
   *  Mesh location structure initialization
   */
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create(1,
                                                           comm,
                                                           PDM_OWNERSHIP_KEEP);

  /* Set target point cloud */
  PDM_mesh_location_n_part_cloud_set(mesh_loc,
                                     0,
                                     1);

  PDM_mesh_location_cloud_set(mesh_loc,
                              0,
                              0,
                              n_pts,
                              pts_coord,
                              pts_ln_to_gn);

  /* Set source mesh */
  PDM_mesh_location_shared_nodal_mesh_set(mesh_loc,
                                          src_pmn);

  /* Set location parameters */
  PDM_mesh_location_tolerance_set(mesh_loc,
                                  tolerance);

  PDM_mesh_location_method_set(mesh_loc,
                               loc_method);

  /*
   *  Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }
  PDM_mesh_location_compute(mesh_loc);

  PDM_mesh_location_dump_times(mesh_loc);

  /*
   *  Check interpolation
   */
  double **src_field = malloc(sizeof(double *) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(src_pmn, ipart);
    double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(src_pmn, ipart);

    src_field[ipart] = malloc(sizeof(double) * n_vtx);
    for (int i = 0; i < n_vtx; i++) {
      src_field[ipart][i] = _eval_field(vtx_coord[3*i  ],
                                        vtx_coord[3*i+1],
                                        vtx_coord[3*i+2],
                                        order_field);
    }
  }

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
                                        0,
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

    int n_elt = PDM_part_mesh_nodal_n_elmts_get(src_pmn,
                                                PDM_GEOMETRY_KIND_VOLUMIC,
                                                ipart);

    send_field[ipart] = malloc(sizeof(double) * elt_pts_idx[n_elt]);
    for (int ielt = 0; ielt < n_elt; ielt++) {
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


  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);

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
  PDM_part_to_part_free(ptp);

  double *tgt_field_interp = malloc(sizeof(double) * n_pts);
  double *tgt_field_exact  = malloc(sizeof(double) * n_pts);
  for (int i = 0; i < n_pts; i++) {
    tgt_field_interp[i] = 123456789;
    tgt_field_exact[i] = _eval_field(pts_coord[3*i  ],
                                     pts_coord[3*i+1],
                                     pts_coord[3*i+2],
                                     order_field);
  }

  int n_located = PDM_mesh_location_n_located_get(mesh_loc,
                                                  0,
                                                  0);

  int *located = PDM_mesh_location_located_get(mesh_loc,
                                               0,
                                               0);

  PDM_g_num_t *location;
  double      *dist2;
  double      *projected_coord;
  PDM_mesh_location_point_location_get(mesh_loc,
                                       0,
                                       0,
                                       &location,
                                       &dist2,
                                       &projected_coord);

  double err_max = 0.;
  for (int i = 0; i < n_located; i++) {
    int pt_id = located[i] - 1;
    tgt_field_interp[pt_id] = recv_field[0][i];

    double err = PDM_ABS(tgt_field_interp[pt_id] - tgt_field_exact[pt_id]);

    if (err > 1e-6) {
      printf("error = %e "PDM_FMT_G_NUM" (%f %f %f) located in "PDM_FMT_G_NUM" at distance %e\n",
             err,
             pts_ln_to_gn[pt_id],
             pts_coord[3*pt_id+0],
             pts_coord[3*pt_id+1],
             pts_coord[3*pt_id+2],
             location[pt_id],
             sqrt(PDM_ABS(dist2[pt_id])));
    }

    err_max = PDM_MAX(err_max, err);
  }
  free(recv_field[0]);
  free(recv_field);


  PDM_g_num_t g_n_wrong = 0;
  double gerr_max;
  PDM_MPI_Allreduce(&err_max, &gerr_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  if (i_rank == 0) {
    printf("gerr_max = %e\n", gerr_max);
    fflush(stdout);
  }


  if (post) {
    char filename[999];
    sprintf(filename, "mesh_location_ho_pts_%d.vtk", i_rank);

    const char   *field_name [2] = {"field_exact", "field_interp"};
    const double *field_value[2] = {tgt_field_exact, tgt_field_interp};
    PDM_vtk_write_std_elements_double(filename,
                                      n_pts,
                                      pts_coord,
                                      pts_ln_to_gn,
                                      PDM_MESH_NODAL_POINT,
                                      n_pts,
                                      NULL,
                                      pts_ln_to_gn,
                                      2,
                                      field_name,
                                      field_value);

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
    for (int ipart = 0; ipart < n_part; ipart++) {
      sprintf(filename, "mesh_location_ho_mesh_part%d_%d.vtk", ipart, i_rank);
      int n_elt = PDM_part_mesh_nodal_section_n_elt_get(src_pmn, 0, ipart);

      int n_vtx = PDM_part_mesh_nodal_n_vtx_get(src_pmn, ipart);
      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(src_pmn, ipart);
      PDM_g_num_t *vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(src_pmn, ipart);

      int         *connec;
      PDM_g_num_t *numabs;
      int         *parent_num;
      PDM_g_num_t *parent_entity_g_num;
      int          _order;
      const char  *ho_ordering;
      PDM_part_mesh_nodal_section_std_ho_get(src_pmn,
                                             0,
                                             ipart,
                                             &connec,
                                             &numabs,
                                             &parent_num,
                                             &parent_entity_g_num,
                                             &_order,
                                             &ho_ordering,
                                             PDM_OWNERSHIP_KEEP);

      int *pcell_vtx_out = malloc(n_vtx_per_elmt * n_elt * sizeof(int));
      for(int i = 0; i < n_vtx_per_elmt * n_elt; ++i) {
        pcell_vtx_out[i] = connec[i];
      }

      if (PDM_Mesh_nodal_elmt_is_ho(elt_type)) {
        PDM_Mesh_nodal_reorder_elt_vtx(elt_type,
                                       _order,
                                       ho_ordering,
                                       "PDM_HO_ORDERING_VTK",
                                       n_elt,
                                       connec,
                                       pcell_vtx_out);
      }

      field_value[0] = src_field[ipart];
      PDM_vtk_write_std_elements_ho_with_vtx_field(filename,
                                                   _order,
                                                   n_vtx,
                                                   vtx_coord,
                                                   vtx_ln_to_gn,
                                                   elt_type,
                                                   n_elt,
                                                   pcell_vtx_out,
                                                   numabs,
                                                   0,
                                                   NULL,
                                                   NULL,
                                                   1,
                                                   field_name,
                                                   field_value);
      free(pcell_vtx_out);
    }
  }

  /* Free memory */
  PDM_mesh_location_free(mesh_loc);
  PDM_part_mesh_nodal_free(src_pmn);
  PDM_multipart_free(src_mpart);
  free(pts_coord);
  free(pts_ln_to_gn);

  free(tgt_field_interp);
  free(tgt_field_exact);

  for (int i = 0; i < n_part; i++) {
    free(send_field[i]);
    free(src_field [i]);
  }
  free(send_field);
  free(src_field );

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return g_n_wrong;
}

