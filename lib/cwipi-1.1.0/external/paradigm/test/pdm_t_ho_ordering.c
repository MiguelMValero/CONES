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
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_mesh_nodal.h"
#include "pdm_ho_ordering.h"
#include "pdm_sort.h"

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
     "  -o     Order\n\n"
     "  -h     This message.\n\n");

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
           int           *order,
           int           *t_elt)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

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

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



// static void _check_ijk_to_ijk
// (
//  const PDM_Mesh_nodal_elt_t  t_elt,
//  const int                   order,
//  int                        *ijk_to_user
//  )
// {
//   int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
//   for (int i = 0; i < n_nodes; i++) {
//     //assert(ijk_to_user[i] == i);
//     if (ijk_to_user[i] != i) {
//       printf("error for t_elt %d\n", (int) t_elt);
//       return;
//     }
//   }
// }

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init(&argc, &argv);

  PDM_ho_ordering_init();

  /*
   *  Set default values
   */
  int order = 1;
  PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_TRIAHO;
  //  1 -> bar
  //  2 -> tria
  //  3 -> quad
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  // 10 -> bar_ho
  // 11 -> tria_ho
  // 12 -> quad_ho
  // 13 -> tetra_ho
  // 14 -> pyramid_ho
  // 15 -> prism_ho
  // 16 -> hexa_ho

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &order,
             (int *) &t_elt);


  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);

  int elt_dim = PDM_Mesh_nodal_elt_dim_get(t_elt);

  int *ijk_to_ijk  = malloc (sizeof(int) * n_nodes * elt_dim);
  int *user_to_ijk = malloc (sizeof(int) * n_nodes * elt_dim);
  int idx = 0;

  switch(t_elt) {
  case PDM_MESH_NODAL_BAR2:
  case PDM_MESH_NODAL_BARHO:
    for (int i = 0; i <= order; i++) {
      ijk_to_ijk[idx++] = i;
    }
    break;

  case PDM_MESH_NODAL_TRIA3:
  case PDM_MESH_NODAL_TRIAHO:
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order-j; i++) {
        ijk_to_ijk[idx++] = i;
        ijk_to_ijk[idx++] = j;
      }
    }
    break;

  case PDM_MESH_NODAL_QUAD4:
  case PDM_MESH_NODAL_QUADHO:
    for (int j = 0; j <= order; j++) {
      for (int i = 0; i <= order; i++) {
        ijk_to_ijk[idx++] = i;
        ijk_to_ijk[idx++] = j;
      }
    }
    break;

  case PDM_MESH_NODAL_TETRA4:
  case PDM_MESH_NODAL_TETRAHO:
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        for (int i = 0; i <= order-j-k; i++) {
          ijk_to_ijk[idx++] = i;
          ijk_to_ijk[idx++] = j;
          ijk_to_ijk[idx++] = k;
        }
      }
    }
    break;

  case PDM_MESH_NODAL_PYRAMID5:
  case PDM_MESH_NODAL_PYRAMIDHO:
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order-k; j++) {
        for (int i = 0; i <= order-k; i++) {
          ijk_to_ijk[idx++] = i;
          ijk_to_ijk[idx++] = j;
          ijk_to_ijk[idx++] = k;
        }
      }
    }
    break;

  case PDM_MESH_NODAL_PRISM6:
  case PDM_MESH_NODAL_PRISMHO:
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        for (int i = 0; i <= order-j; i++) {
          ijk_to_ijk[idx++] = i;
          ijk_to_ijk[idx++] = j;
          ijk_to_ijk[idx++] = k;
        }
      }
    }
    break;

  case PDM_MESH_NODAL_HEXA8:
  case PDM_MESH_NODAL_HEXAHO:
    for (int k = 0; k <= order; k++) {
      for (int j = 0; j <= order; j++) {
        for (int i = 0; i <= order; i++) {
          ijk_to_ijk[idx++] = i;
          ijk_to_ijk[idx++] = j;
          ijk_to_ijk[idx++] = k;
        }
      }
    }
    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt %d\n", (int) t_elt);
  }

  int n_rep = 5;
  char name[999];


  double *rnd  = malloc (sizeof(double) * n_nodes);
  int    *perm = malloc (sizeof(int   ) * n_nodes);

  for (int irep = 0; irep < n_rep; irep++) {

    for (int i = 0; i < n_nodes; i++) {
      perm[i] = i;
      rnd[i] = rand();
    }

    PDM_sort_double (rnd, perm, n_nodes);

    for (int i = 0; i < n_nodes; i++) {
      int inode = perm[i];

      for (int j = 0; j < elt_dim; j++) {
        user_to_ijk[elt_dim*i + j] = ijk_to_ijk[elt_dim*inode + j];
      }
    }

    sprintf(name, "MY_ORDERING_%d", irep);

    PDM_ho_ordering_user_to_ijk_add (name,
                                     t_elt,
                                     order,
                                     n_nodes,
                                     user_to_ijk);
  }

  free (perm);
  free (rnd);
  free (ijk_to_ijk);
  free (user_to_ijk);


  /*int *user_to_ijk = PDM_ho_ordering_user_to_ijk_get ("PDM_HO_ORDERING_VTK",
                                                      t_elt,
                                                      order);
  int *ijk_to_user = PDM_ho_ordering_ijk_to_user_get ("PDM_HO_ORDERING_VTK",
                                                      t_elt,
                                                      order);

  if (user_to_ijk != NULL) {
    for (int i = 0; i < n_nodes; i++) {
      printf("user node %2d: ", i);
      for (int j = 0; j < elt_dim; j++) {
        printf("%2d ", user_to_ijk[elt_dim*i + j]);
      }
      printf("\n");
    }
  }

  if (ijk_to_user != NULL) {
    for (int i = 0; i < n_nodes; i++) {
      printf("ijk node %3d --> user node %3d\n", i, ijk_to_user[i]);
    }
    }*/

  /*int n_nodes_max = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_HEXA8, order);

  int *user_to_ijk = malloc (sizeof(int) * n_nodes_max * 3);
  int *ijk_to_user = NULL;
  int idx;


  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_BAR2;
       type <= PDM_MESH_NODAL_HEXA8; type++) {

    idx = 0;

    switch(type) {
    case PDM_MESH_NODAL_BAR2:
      for (int i = 0; i <= order; i++) {
        user_to_ijk[idx++] = i;
      }
      break;

    case PDM_MESH_NODAL_TRIA3:
      for (int j = 0; j <= order; j++) {
        for (int i = 0; i <= order-j; i++) {
          user_to_ijk[idx++] = i;
          user_to_ijk[idx++] = j;
        }
      }
      break;

    case PDM_MESH_NODAL_QUAD4:
      for (int j = 0; j <= order; j++) {
        for (int i = 0; i <= order; i++) {
          user_to_ijk[idx++] = i;
          user_to_ijk[idx++] = j;
        }
      }
      break;

      case PDM_MESH_NODAL_TETRA4:
      for (int k = 0; k <= order; k++) {
        for (int j = 0; j <= order-k; j++) {
          for (int i = 0; i <= order-j-k; i++) {
            user_to_ijk[idx++] = i;
            user_to_ijk[idx++] = j;
            user_to_ijk[idx++] = k;
          }
        }
      }
      break;

    case PDM_MESH_NODAL_PYRAMID5:
      if (order < 2) {
        for (int k = 0; k <= order; k++) {
          for (int j = 0; j <= order-k; j++) {
            for (int i = 0; i <= order-k; i++) {
              user_to_ijk[idx++] = i;
              user_to_ijk[idx++] = j;
              user_to_ijk[idx++] = k;
            }
          }
        }
      } else {
        continue;
      }
      break;

    case PDM_MESH_NODAL_PRISM6:
      for (int k = 0; k <= order; k++) {
        for (int j = 0; j <= order; j++) {
          for (int i = 0; i <= order-j; i++) {
            user_to_ijk[idx++] = i;
            user_to_ijk[idx++] = j;
            user_to_ijk[idx++] = k;
          }
        }
      }
      break;

    case PDM_MESH_NODAL_HEXA8:
      for (int k = 0; k <= order; k++) {
        for (int j = 0; j <= order; j++) {
          for (int i = 0; i <= order; i++) {
            user_to_ijk[idx++] = i;
            user_to_ijk[idx++] = j;
            user_to_ijk[idx++] = k;
          }
        }
      }
      break;

    default:
      continue;
    }

    ijk_to_user = PDM_ho_ordering_compute_ijk_to_user (type,
                                                       order,
                                                       user_to_ijk);

    log_trace("type = %d\n", (int) (type));
    PDM_log_trace_array_int (ijk_to_user,
                             PDM_Mesh_nodal_n_vtx_elt_get (type, order),
                             "ijk_to_user : ");

    _check_ijk_to_ijk (type, order, ijk_to_user);

    printf("type %d OK\n", (int) (type));

    free (ijk_to_user);
  }

  free (user_to_ijk);*/


  PDM_MPI_Finalize();

  return 0;
}
