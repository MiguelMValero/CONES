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
#include "pdm_dcube_nodal_gen.h"
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
  PDM_generate_lines(comm,
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
  PDM_dmesh_t* dm = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                     0,
                                     0,
                                     dn_edge,
                                     dn_vtx,
                                     comm);

  PDM_dmesh_connectivity_set(dm,
                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                             dedge_vtx,
                             NULL,
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_vtx_coord_set(dm,
                          dvtx_coord,
                          PDM_OWNERSHIP_USER);

  /*
   * Mulitpart
   */
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  // PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
  int n_part = 1;
  PDM_multipart_t* mpart = PDM_multipart_create(1,
                                                &n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_dmesh_set(mpart, 0, dm);

  PDM_multipart_compute(mpart);

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
  PDM_dmesh_free(dm);

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
