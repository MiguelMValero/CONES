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
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_octree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_multipart.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

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
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *n_g_pts
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_g_pts = atol(argv[i]);
        *n_g_pts = (PDM_g_num_t) _n_g_pts;
      }
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

static
void
_generate_lines
(
  PDM_MPI_Comm  comm,
  double        zero_x,
  double        zero_y,
  double        zero_z,
  double        length,
  PDM_g_num_t   n_g_pts,
  PDM_g_num_t **distrib_edge_out,
  PDM_g_num_t **distrib_vtx_out,
  PDM_g_num_t **dedge_vtx_out,
  double      **dvtx_coord_out
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t gn_vtx  = (n_g_pts    );
  PDM_g_num_t gn_edge = (n_g_pts - 1);

  int dcube_nx = n_g_pts - 1;

  PDM_g_num_t* distrib_edge = PDM_compute_uniform_entity_distribution(comm, gn_edge);
  PDM_g_num_t* distrib_vtx  = PDM_compute_uniform_entity_distribution(comm, gn_vtx);

  int dn_vtx  = (int) (distrib_vtx [i_rank+1] - distrib_vtx [i_rank]);
  int dn_edge = (int) (distrib_edge[i_rank+1] - distrib_edge[i_rank]);

  double *dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);

  double step_x = length / (double) (n_g_pts - 1);

  for (int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {

    PDM_g_num_t g_vtx = distrib_vtx[i_rank] + i_vtx;

    PDM_g_num_t indi = g_vtx % n_g_pts;

    dvtx_coord[3 * i_vtx    ] = indi * step_x + zero_x;
    dvtx_coord[3 * i_vtx + 1] = zero_y;
    dvtx_coord[3 * i_vtx + 2] = zero_z;
  }


  PDM_g_num_t *dedge_vtx     = malloc( 2 * dn_edge * sizeof(PDM_g_num_t));

  for (int i_edge = 0; i_edge < dn_edge; ++i_edge) {

    PDM_g_num_t g = distrib_edge[i_rank] + i_edge;

    PDM_g_num_t indi = g % dcube_nx;

    dedge_vtx[2*i_edge  ] = 1 + (indi  );
    dedge_vtx[2*i_edge+1] = 1 + (indi+1);
  }

  *dvtx_coord_out   = dvtx_coord;
  *dedge_vtx_out    = dedge_vtx;
  *distrib_edge_out = distrib_edge;
  *distrib_vtx_out  = distrib_vtx;
}




/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t n_g_pts   = 10;
  _read_args(argc,
             argv,
             &n_g_pts);

  double      *dvtx_coord   = NULL;
  PDM_g_num_t *distrib_edge = NULL;
  PDM_g_num_t *distrib_vtx  = NULL;
  PDM_g_num_t *dedge_vtx    = NULL;
  _generate_lines(comm,
                  0.,
                  0.,
                  0.,
                  1.,
                  n_g_pts,
                  &distrib_edge,
                  &distrib_vtx,
                  &dedge_vtx,
                  &dvtx_coord);

  int dn_edge = distrib_edge[i_rank+1] - distrib_edge[i_rank];
  int dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank];

  /*
   * Create dmesh
   */
  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(comm,
                                                  1,
                                                  distrib_vtx [n_rank],
                                                  0,
                                                  0,
                                                  distrib_edge[n_rank]);

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_USER);

  int id_section_edge = PDM_DMesh_nodal_section_add(dmn,
                                                    PDM_GEOMETRY_KIND_RIDGE,
                                                    PDM_MESH_NODAL_BAR2);
  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_RIDGE,
                                  id_section_edge,
                                  dn_edge,
                                  dedge_vtx,
                                  PDM_OWNERSHIP_USER);

  /*
   * Mulitpart
   */
  // PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
  int n_part = 1;
  PDM_multipart_t* mpart = PDM_multipart_create(1,
                                                &n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

  PDM_multipart_compute(mpart);
  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmn, PDM_OWNERSHIP_KEEP);
  PDM_part_mesh_nodal_free(pmn);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      PDM_g_num_t* pvtx_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &pvtx_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* pedge_ln_to_gn = NULL;
      int pn_edge = PDM_multipart_part_ln_to_gn_get(mpart,
                                                    0,
                                                    i_part,
                                                    PDM_MESH_ENTITY_EDGE,
                                                    &pedge_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);

      int *pedge_vtx     = NULL;
      int *pedge_vtx_idx = NULL;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &pedge_vtx_idx,
                                          &pedge_vtx,
                                          PDM_OWNERSHIP_KEEP);

      double *pvtx_coord = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &pvtx_coord,
                                                   PDM_OWNERSHIP_KEEP);
      char filename[999];
      sprintf(filename, "out_part_vtx_%i.vtk", i_rank);
      PDM_vtk_write_std_elements (filename,
                                  n_vtx,
                                  pvtx_coord,
                                  pvtx_ln_to_gn,
                                  PDM_MESH_NODAL_BAR2,
                                  pn_edge,
                                  pedge_vtx,
                                  pedge_ln_to_gn,
                                  0,
                                  NULL,
                                  NULL);

    }
  }


  PDM_multipart_free(mpart);
  PDM_DMesh_nodal_free(dmn);

  free (dvtx_coord);
  free (distrib_vtx);
  free (distrib_edge);
  free (dedge_vtx);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
