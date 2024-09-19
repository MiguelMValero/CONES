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
#include "pdm_dist_cloud_surf.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_multipart.h"
#include "pdm_block_to_block.h"
#include "pdm_part_to_block.h"
#include "pdm_distrib.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_reader_stl.h"

#include "pdm_vtk.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
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
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_vtx_seg,
           char         **filename,
           PDM_g_num_t   *n_g_pts_clouds,
           double        *length,
           int           *n_part,
           int           *post,
           int           *method)
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
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_g_pts_clouds = atol(argv[i]);
        *n_g_pts_clouds = (PDM_g_num_t) _n_g_pts_clouds;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = PDM_PART_SPLIT_HILBERT;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static
void
_split_surface_mesh
(
 const PDM_MPI_Comm        comm,
 const PDM_split_dual_t    part_method,
 const int                 n_part,
       PDM_dmesh_nodal_t  *dmn,
       PDM_multipart_t   **_mpart
)
{
  int n_domain = 1;
  int *n_part_domains = (int *) malloc(sizeof(int) * n_domain);
  n_part_domains[0] = n_part;

  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                n_part_domains,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);

  free(n_part_domains);

  *_mpart = mpart;
}



int main(int argc, char *argv[])
{
  /*
   *  Init
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Set default values
   */

  PDM_g_num_t           n_vtx_seg      = 100;
  PDM_g_num_t           n_g_pts_clouds = 10;
  double                length         = 1.;
  int                   n_part         = 1;
  int                   post           = 0;
  char                 *filename  = NULL;

  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_HILBERT;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &filename,
             &n_g_pts_clouds,
             &length,
             &n_part,
             &post,
     (int *) &part_method);

  if (filename == NULL) {
    filename = (char *) PDM_MESH_DIR"sphere.stl";
  }

  double radius = length;//2*length;

  /*
   *  Generate cloud
   */
  if (i_rank == 0) {
    printf("-- Generate volume mesh\n");
    fflush(stdout);
  }
  int          n_pts_clouds;
  double      *pts_coord;
  PDM_g_num_t *pts_g_num;
  PDM_point_cloud_gen_random (comm,
                              0, // seed
                              0, // geometric_g_num
                              n_g_pts_clouds,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_pts_clouds,
                              &pts_coord,
                              &pts_g_num);


  /*
   *  Read surface mesh
   */
  if (i_rank == 0) {
    printf("-- Read surface mesh : %s \n", filename);
    fflush(stdout);
  }

  PDM_dmesh_nodal_t* dmn = PDM_reader_stl_dmesh_nodal(comm,
                                                      filename);

  PDM_multipart_t *mpart_surf = NULL;
  _split_surface_mesh(comm,
                      part_method,
                      n_part,
                      dmn,
                      &mpart_surf);

  int          *surf_pn_vtx          = (int          *) malloc(sizeof(int          ) * n_part);
  int          *surf_pn_face         = (int          *) malloc(sizeof(int          ) * n_part);
  int          *surf_pn_edge         = (int          *) malloc(sizeof(int          ) * n_part);
  int         **surf_pface_edge_idx  = (int         **) malloc(sizeof(int         *) * n_part);
  int         **surf_pface_edge      = (int         **) malloc(sizeof(int         *) * n_part);
  int         **surf_pedge_vtx       = (int         **) malloc(sizeof(int         *) * n_part);
  int         **surf_pface_vtx       = (int         **) malloc(sizeof(int         *) * n_part);
  double      **surf_pvtx_coord      = (double      **) malloc(sizeof(double      *) * n_part);
  PDM_g_num_t **surf_pvtx_ln_to_gn   = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **surf_pface_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    surf_pn_vtx[i_part] = PDM_multipart_part_vtx_coord_get(mpart_surf,
                                                           0,
                                                           i_part,
                                                           &surf_pvtx_coord[i_part],
                                                           PDM_OWNERSHIP_KEEP);

    PDM_multipart_part_ln_to_gn_get(mpart_surf,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    &surf_pvtx_ln_to_gn[i_part],
                                    PDM_OWNERSHIP_KEEP);

    surf_pn_face[i_part] = PDM_multipart_part_ln_to_gn_get(mpart_surf,
                                                           0,
                                                           i_part,
                                                           PDM_MESH_ENTITY_FACE,
                                                           &surf_pface_ln_to_gn[i_part],
                                                           PDM_OWNERSHIP_KEEP);

    surf_pn_face[i_part] = PDM_multipart_part_connectivity_get(mpart_surf,
                                                               0,
                                                               i_part,
                                                               PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                               &surf_pface_edge_idx[i_part],
                                                               &surf_pface_edge[i_part],
                                                               PDM_OWNERSHIP_KEEP);

    int* tmp_pedge_vtx_idx = NULL;
    surf_pn_edge[i_part] = PDM_multipart_part_connectivity_get(mpart_surf,
                                                               0,
                                                               i_part,
                                                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                               &tmp_pedge_vtx_idx,
                                                               &surf_pedge_vtx[i_part],
                                                               PDM_OWNERSHIP_KEEP);
    assert(tmp_pedge_vtx_idx == NULL);


    /*
     * Compute face_vtx
     */
    PDM_compute_face_vtx_from_face_and_edge(surf_pn_face       [i_part],
                                            surf_pface_edge_idx[i_part],
                                            surf_pface_edge    [i_part],
                                            surf_pedge_vtx     [i_part],
                                            &surf_pface_vtx    [i_part]);

    if(0 == 1) {
      PDM_log_trace_array_long(surf_pface_ln_to_gn[i_part], surf_pn_face[i_part]                            , "surf_pface_ln_to_gn :: ");
      PDM_log_trace_array_int (surf_pface_edge_idx[i_part], surf_pn_face[i_part]+1                          , "surf_pface_edge_idx :: ");
      // PDM_log_trace_array_int(surf_pface_vtx[i_part]     , surf_pface_edge_idx[i_part][surf_pn_face[i_part]], "surf_pface_vtx     :: ");
    }
  }


  /*
   * Identity for each point in cloud if it's inside or outside surf
   */
  PDM_dist_cloud_surf_t* ics = PDM_dist_cloud_surf_create(PDM_MESH_NATURE_MESH_SETTED,
                                                          1,
                                                          comm,
                                                          PDM_OWNERSHIP_KEEP);

  int n_part_cloud  = 1;
  int i_point_cloud = 0;
  int i_part_cloud  = 0;
  PDM_dist_cloud_surf_n_part_cloud_set(ics, i_part_cloud, n_part_cloud);

  PDM_dist_cloud_surf_cloud_set(ics,
                                  i_point_cloud,
                                  i_part_cloud,
                                  n_pts_clouds,
                                  pts_coord,
                                  pts_g_num);

  PDM_dist_cloud_surf_surf_mesh_global_data_set(ics, n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_dist_cloud_surf_surf_mesh_part_set(ics,
                                             i_part,
                                             surf_pn_face       [i_part],
                                             surf_pface_edge_idx[i_part],
                                             surf_pface_vtx     [i_part],
                                             surf_pface_ln_to_gn[i_part],
                                             surf_pn_vtx        [i_part],
                                             surf_pvtx_coord    [i_part],
                                             surf_pvtx_ln_to_gn [i_part]);
  }

  PDM_dist_cloud_surf_compute(ics);

  PDM_dist_cloud_surf_dump_times(ics);

  double      *distance = NULL;
  double      *projected = NULL;
  PDM_g_num_t *closest_elt_gnum = NULL;
  PDM_dist_cloud_surf_get(ics,
                          i_point_cloud,
                          i_part_cloud,
                          &distance,
                          &projected,
                          &closest_elt_gnum);


  // if (1 == 1) {
  //   for (int i_part = 0; i_part < n_part_cloud; ++i_part)
  //   {
  //     char filename_out[999];

  //     sprintf(filename_out, "raytracing_distance_%2.2d.vtk", i_rank);
  //     PDM_vtk_write_point_cloud(filename_out,
  //                               n_pts_clouds,
  //                               pts_coord,
  //                               NULL,
  //                               distance);
  //                               // &distance[0],
  //                               // NULL);
  //   }
  // }

  PDM_dist_cloud_surf_free(ics);

  free(pts_coord);
  free(pts_g_num);
  free(surf_pn_edge       );
  free(surf_pn_vtx        );
  free(surf_pn_face       );
  free(surf_pface_edge    );
  free(surf_pface_edge_idx);
  for(int i_part = 0; i_part < n_part; i_part++) {
    free(surf_pface_vtx[i_part]);
  }
  free(surf_pface_vtx     );
  free(surf_pedge_vtx     );
  free(surf_pvtx_coord    );
  free(surf_pvtx_ln_to_gn );
  free(surf_pface_ln_to_gn);


  PDM_DMesh_nodal_free(dmn);

  PDM_multipart_free(mpart_surf);

  PDM_MPI_Barrier (comm);
  if (i_rank == 0) {
  PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize ();

  return 0;
}
