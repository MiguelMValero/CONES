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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_geom_elem.h"
#include "pdm_priv.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_ho_location.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_array.h"
#include "pdm_ho_basis.h"


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
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
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
_read_args
(
 int                    argc,
 char                 **argv,
 int                   *n_pts,
 int                   *order,
 PDM_Mesh_nodal_elt_t  *t_elt,
 int                   *visu,
 int                   *seed
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_pts = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-o") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *order = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }
    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *seed = atoi(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static int
_gen_ho_elt
(
 const PDM_Mesh_nodal_elt_t   t_elt,
 const int                    order,
       double               **node_coord
)
{
  int n_node = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

  *node_coord = malloc(sizeof(double) * n_node * 3);

  switch (t_elt) {
    case PDM_MESH_NODAL_BARHO:
    case PDM_MESH_NODAL_BARHO_BEZIER: {
      for (int i = 0; i < n_node; i++) {
        double x = (double) i / (double) order - 0.5;
        (*node_coord)[3*i  ] = x;
        (*node_coord)[3*i+1] = 0.2*x*x;
        (*node_coord)[3*i+2] = 0.1*x*x;
      }
      break;
    }

    case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER: {
      int idx = 0;
      for (int j = 0; j <= order; j++) {
        double y = (double) j / (double) order - 0.5;
        for (int i = 0; i <= order - j; i++) {
          double x = (double) i / (double) order - 0.5;
          (*node_coord)[idx++] = x;
          (*node_coord)[idx++] = y;
          (*node_coord)[idx++] = 0.1*(x*x + y*y);
        }
      }
      break;
    }

    case PDM_MESH_NODAL_QUADHO: {
      int idx = 0;
      for (int j = 0; j <= order; j++) {
        double y = (double) j / (double) order - 0.5;
        for (int i = 0; i <= order; i++) {
          double x = (double) i / (double) order - 0.5;
          (*node_coord)[idx++] = x;
          (*node_coord)[idx++] = y;
          (*node_coord)[idx++] = 0.1*(x*x + y*y);
        }
      }
      break;
    }

    case PDM_MESH_NODAL_TETRAHO: {
      int idx = 0;
      for (int k = 0; k <= order; k++) {
        double z = (double) k / (double) order - 0.5;
        for (int j = 0; j <= order-k; j++) {
          double y = (double) j / (double) order - 0.5;
          for (int i = 0; i <= order-j-k; i++) {
            double x = (double) i / (double) order - 0.5;
            (*node_coord)[idx++] = x + 0.2*z*z;
            (*node_coord)[idx++] = y - 0.2*x*x;
            (*node_coord)[idx++] = z + 0.2*y*y;
          }
        }
      }
      break;
    }

    case PDM_MESH_NODAL_PYRAMIDHO: {
      int idx = 0;
      for (int k = 0; k <= order; k++) {
        double z = (double) k / (double) order - 0.5;
        for (int j = 0; j <= order-k; j++) {
          double y = (double) j / (double) order - 0.5;
          for (int i = 0; i <= order-k; i++) {
            double x = (double) i / (double) order - 0.5;
            (*node_coord)[idx++] = x + 0.2*z*z;
            (*node_coord)[idx++] = y - 0.2*x*x;
            (*node_coord)[idx++] = z + 0.2*y*y;
          }
        }
      }
      break;
    }

    case PDM_MESH_NODAL_PRISMHO: {
      int idx = 0;
      for (int k = 0; k <= order; k++) {
        double z = (double) k / (double) order - 0.5;
        for (int j = 0; j <= order; j++) {
          double y = (double) j / (double) order - 0.5;
          for (int i = 0; i <= order-j; i++) {
            double x = (double) i / (double) order - 0.5;
            (*node_coord)[idx++] = x + 0.2*z*z;
            (*node_coord)[idx++] = y - 0.2*x*x;
            (*node_coord)[idx++] = z + 0.2*y*y;
          }
        }
      }
      break;
    }

    case PDM_MESH_NODAL_HEXAHO: {
      int idx = 0;
      for (int k = 0; k <= order; k++) {
        double z = (double) k / (double) order - 0.5;
        for (int j = 0; j <= order; j++) {
          double y = (double) j / (double) order - 0.5;
          for (int i = 0; i <= order; i++) {
            double x = (double) i / (double) order - 0.5;
            (*node_coord)[idx++] = x + 0.2*z*z;
            (*node_coord)[idx++] = y - 0.2*x*x;
            (*node_coord)[idx++] = z + 0.2*y*y;
          }
        }
      }
      break;
    }

    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid elt_type %d\n", (int) t_elt);

  }

  return n_node;
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
  int                  n_pts = 1;
  int                  order = 3;
  int                  visu  = 0;
  PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_QUADHO;
  int                  seed  = -1;
  // 11 -> tria_ho
  // 12 -> quad_ho
  // 13 -> tetra_ho
  // 14 -> pyramid_ho
  // 15 -> prism_ho
  // 16 -> hexa_ho

  assert(t_elt > PDM_MESH_NODAL_POLY_3D);


  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_pts,
             &order,
             &t_elt,
             &visu,
             &seed);

  PDM_MPI_Init(&argc, &argv);


  double *parent_node_coord = NULL;
  int n_node = _gen_ho_elt(t_elt,
                           order,
                           &parent_node_coord);


  int elt_dim = PDM_Mesh_nodal_elt_dim_get(t_elt);




  double uvw[3];

  // uvw[0] = -0.1;
  // uvw[1] = 1.2;
  // uvw[2] = 1.3;
  if (seed < 0) {
    seed = time(NULL);
  }
  printf("seed = %d\n", seed);
  srand(seed);
  for (int i = 0; i < elt_dim; i++) {
    uvw[i] = 0.5 + 1.0*(2*((double) rand() / (double) RAND_MAX) - 1);
  }

  double *pts_coord = malloc(sizeof(double) * 3 * n_pts);
  double *weight = malloc(sizeof(double) * n_node);
  PDM_ho_basis(t_elt,
               order,
               n_node,
               1,
               uvw,
               weight);

  for (int j = 0; j < 3; j++) {
    pts_coord[j] = 0;
  }

  for (int i = 0; i < n_node; i++) {
    for (int j = 0; j < 3; j++) {
      pts_coord[j] += weight[i] * parent_node_coord[3*i+j];
    }
  }

  const double tolerance = 1e-6;
  double *work_array = malloc(sizeof(double) * n_node * 4);
  double uvw2[3];
  double dist2 = 0;
  double *proj_coord = malloc(sizeof(double) * n_pts * 3);
  int converged;
  dist2 = PDM_ho_location_newton(t_elt,
                                 order,
                                 n_node,
                                 parent_node_coord,
                                 pts_coord,
                                 tolerance,
                                 proj_coord,
                                 uvw2,
                                 &converged,
                                 work_array);
  printf("converged? %d\n", converged);
  free(work_array);


  double uvw3[3];
  double *proj_coord3 = malloc(sizeof(double) * n_pts * 3);
  double dist3 = PDM_ho_location(t_elt,
                                 order,
                                 n_node,
                                 parent_node_coord,
                                 pts_coord,
                                 proj_coord3,
                                 uvw3);

  printf("uvw  : ");
  for (int i = 0; i < elt_dim; i++) {
    printf("%e ", uvw[i]);
  }
  printf("\n");

  printf("uvw2 : ");
  for (int i = 0; i < elt_dim; i++) {
    printf("%e ", uvw2[i]);
  }
  printf("\n");

  printf("uvw3 : ");
  for (int i = 0; i < elt_dim; i++) {
    printf("%e ", uvw3[i]);
  }
  printf("\n");





  printf("dist_newton = %e, dist_subdiv = %e, delta = %e, relatif = %e\n",
         sqrt(dist2), sqrt(dist3), sqrt(dist2) - sqrt(dist3), (sqrt(dist2) - sqrt(dist3))/sqrt(dist2));

  free(weight);


  if (visu) {
    PDM_vtk_write_std_elements("point.vtk",
                               1,
                               pts_coord,
                               NULL,
                               PDM_MESH_NODAL_POINT,
                               1,
                               NULL,
                               NULL,
                               0,
                               NULL,
                               NULL);

    PDM_vtk_write_std_elements("proj.vtk",
                               1,
                               proj_coord,
                               NULL,
                               PDM_MESH_NODAL_POINT,
                               1,
                               NULL,
                               NULL,
                               0,
                               NULL,
                               NULL);

    PDM_vtk_write_std_elements("proj_subdiv.vtk",
                               1,
                               proj_coord3,
                               NULL,
                               PDM_MESH_NODAL_POINT,
                               1,
                               NULL,
                               NULL,
                               0,
                               NULL,
                               NULL);

    int *connec = malloc(sizeof(int) * n_node);

    int *ijk_to_user = PDM_ho_ordering_ijk_to_user_get("PDM_HO_ORDERING_VTK",
                                                       t_elt,
                                                       order);

    for (int i = 0; i < n_node; i++) {
      connec[ijk_to_user[i]] = i + 1;
    }

    PDM_vtk_write_std_elements_ho("parent.vtk",
                                  order,
                                  n_node,
                                  parent_node_coord,
                                  NULL,
                                  t_elt,
                                  1,
                                  connec,
                                  NULL,
                                  0,
                                  NULL,
                                  NULL);

    free(connec);
  }

  PDM_MPI_Finalize();


  free(parent_node_coord);
  free(pts_coord);
  free(proj_coord);
  free(proj_coord3);

  return 0;
}

