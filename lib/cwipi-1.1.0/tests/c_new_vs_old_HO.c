/*
  This file is part of the CWIPI library.

  Copyright (C) 2023  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "cwipi.h"
#include "cwp.h"
#include "cwp_priv.h"

#include "grid_mesh.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dmesh.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_vtk.h"

#define ABS(a)   ((a) <  0  ? -(a) : (a))
#define MAX(a,b) ((a) > (b) ?  (a) : (b))

typedef enum {
  CWP_VERSION_OLD,
  CWP_VERSION_NEW
} CWP_Version_t;

/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -n           <> Number of vertices in band length.\n\n"
         "  -no_random      Disable mesh randomization\n\n"
         "  -n_proc_data <> Number of processes where there are data \n\n"
         "  -h              this message.\n\n");

  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

static void
_read_args
(
  int                    argc,
  char                 **argv,
  CWP_Version_t         *version,
  PDM_g_num_t           *n_vtx_seg,
  int                   *deform,
  double                *tolerance,
  PDM_split_dual_t      *part_method,
  PDM_Mesh_nodal_elt_t  *elt_type,
  int                   *order
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-new") == 0) {
      *version = CWP_VERSION_NEW;
    }
    else if (strcmp(argv[i], "-old") == 0) {
      *version = CWP_VERSION_OLD;
    }
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg = atoi(argv[i]);
      }
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
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
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
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *order = atoi(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

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


static void
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
       int                   **n_node,
       double               ***node_coord,
       PDM_g_num_t          ***node_ln_to_gn,
       int                   **n_elt,
       int                  ***elt_node,
       PDM_g_num_t          ***elt_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_zone = 1;
  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                &n_part,
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

  PDM_dcube_nodal_gen_ordering_set(dcube, "PDM_HO_ORDERING_VTK");

  PDM_dcube_nodal_gen_build(dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  if (deform) {
    int order_deformation = MAX(1, (int) ceil (sqrt(order)));
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
  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);

  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmn, PDM_OWNERSHIP_KEEP);

  *n_node        = malloc(sizeof(int          ) * n_part);
  *n_elt         = malloc(sizeof(int          ) * n_part);
  *node_coord    = malloc(sizeof(double      *) * n_part);
  *elt_node      = malloc(sizeof(int         *) * n_part);
  *node_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);
  *elt_ln_to_gn  = malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    (*n_node)[ipart] = PDM_part_mesh_nodal_n_vtx_get        (pmn, ipart);
    (*n_elt) [ipart] = PDM_part_mesh_nodal_section_n_elt_get(pmn, 0, ipart);

    double *_node_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, ipart);
    (*node_coord)[ipart] = malloc(sizeof(double) * (*n_node)[ipart] * 3);
    memcpy((*node_coord)[ipart], _node_coord, sizeof(double) * (*n_node)[ipart] * 3);

    PDM_g_num_t *_node_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmn, ipart);
    (*node_ln_to_gn)[ipart] = malloc(sizeof(PDM_g_num_t) * (*n_node)[ipart]);
    memcpy((*node_ln_to_gn)[ipart], _node_ln_to_gn, sizeof(PDM_g_num_t) * (*n_node)[ipart]);


    int         *connec;
    PDM_g_num_t *numabs;
    int         *parent_num;
    PDM_g_num_t *parent_entity_g_num;
    int          _order;
    const char  *ho_ordering;
    PDM_part_mesh_nodal_section_std_ho_get(pmn,
                                           0,
                                           ipart,
                                           &connec,
                                           &numabs,
                                           &parent_num,
                                           &parent_entity_g_num,
                                           &_order,
                                           &ho_ordering,
                                           PDM_OWNERSHIP_KEEP);

    int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, _order);

    (*elt_node)[ipart] = malloc(sizeof(int) * n_vtx_elt * (*n_elt)[ipart]);
    memcpy((*elt_node)[ipart], connec, sizeof(int) * n_vtx_elt * (*n_elt)[ipart]);

    (*elt_ln_to_gn)[ipart] = malloc(sizeof(PDM_g_num_t) * (*n_elt)[ipart]);
    memcpy((*elt_ln_to_gn)[ipart], numabs, sizeof(PDM_g_num_t) * (*n_elt)[ipart]);
  }

  PDM_part_mesh_nodal_free(pmn);

  PDM_dcube_nodal_gen_free(dcube);

  PDM_multipart_free(mpart);
}



static int
_ijk_grid
(
 const PDM_Mesh_nodal_elt_t   t_elt,
 const int                    order,
       int                  **ijk
 )
{
  int n_node  = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
  int elt_dim = PDM_Mesh_nodal_elt_dim_get  (t_elt);


  *ijk = malloc(sizeof(int) * n_node * elt_dim);

  switch (t_elt) {
    case PDM_MESH_NODAL_HEXAHO: {
      int idx = 0;
      for (int k = 0; k <= order; k++) {
        for (int j = 0; j <= order; j++) {
          for (int i = 0; i <= order; i++) {
            (*ijk)[idx++] = i;
            (*ijk)[idx++] = j;
            (*ijk)[idx++] = k;
          }
        }
      }
      break;
    }
    default:
      PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt\n");
  }

  return n_node;
}


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

int
main
(
 int   argc,
 char *argv[]
 ) {
  /*
   *  Set default values
   */
  CWP_Version_t        version    = CWP_VERSION_NEW;
  PDM_g_num_t          n_vtx_seg  = 4;
  int                  deform     = 0;
  double               tolerance  = 1e-6;
  int                  n_part     = 1;
  int                  order      = 2;
  PDM_Mesh_nodal_elt_t elt_type   = PDM_MESH_NODAL_TRIAHO; // PDM_MESH_NODAL_HEXAHO;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &version,
             &n_vtx_seg,
             &deform,
             &tolerance,
             &part_method,
             &elt_type,
             &order);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;

  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  assert(n_rank > 1);


  int elt_dim = PDM_Mesh_nodal_elt_dim_get(elt_type);


  // Initialize CWIPI
  const char **code_name         = malloc(sizeof(char *) * 2);
  const char **coupled_code_name = malloc(sizeof(char *) * 2);
  CWP_Status_t is_active_rank    = CWP_STATUS_ON;


  int has_code[2] = {0, 0};
  if (i_rank < n_rank/2) {//(2*n_rank) / 3) {
    has_code[0] = 1;
  }
  if (i_rank >= n_rank/2) {//n_rank / 3) {
    has_code[1] = 1;
  }

  const char *all_code_names[2] = {"code1", "code2"};

  int n_code = 0;
  int code_id  [2];
  for (int icode = 0; icode < 2; icode++) {
    if (has_code[icode]) {
      // log_trace("I run %s\n", all_code_names[icode]);
      code_id          [n_code] = icode+1;
      code_name        [n_code] = all_code_names[icode];
      coupled_code_name[n_code] = all_code_names[(icode+1)%2];
      n_code++;
    }
  }

  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);
  if (version == CWP_VERSION_OLD) {
    assert(n_code == 1);
    assert(n_part == 1);

    cwipi_init(comm,
               code_name[0],
               intra_comm);
  }
  else {
    CWP_Init(comm,
             n_code,
             (const char **) code_name,
             is_active_rank,
             intra_comm);
  }

  if (i_rank == 0) {
    printf("CWIPI Init OK\n");
  }

  // Create coupling
  const char *cpl_name = "c_new_api_vol_HO";
  CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE;
  // CWP_Spatial_interp_t spatial_interp = CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES;

  if (version == CWP_VERSION_OLD) {
    cwipi_create_coupling(cpl_name,                                  // Coupling id
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                          coupled_code_name[0],                      // Coupled application id
                          elt_dim,                                   // Geometric entities dimension
                          tolerance,                                 // Geometric tolerance
                          CWIPI_STATIC_MESH,                         // Mesh type
                          CWIPI_SOLVER_CELL_VERTEX,                  // Solver type
                          -1,                                        // Postprocessing frequency
                          "EnSight Gold",                            // Postprocessing format
                          "text");
  }
  else {
    for (int i_code = 0; i_code < n_code; i_code++) {
      CWP_Cpl_create(code_name[i_code],                                     // Code name
                     cpl_name,                                              // Coupling id
                     coupled_code_name[i_code],                             // Coupled application id
                     (CWP_Interface_t) elt_dim,
                     CWP_COMM_PAR_WITH_PART,                                // Coupling type
                     spatial_interp,
                     n_part,                                                // Number of partitions
                     CWP_DYNAMIC_MESH_STATIC,                               // Mesh displacement type
                     CWP_TIME_EXCH_USER_CONTROLLED);                        // Postprocessing frequency

      // CWP_Visu_set(code_name[i_code],       // Code name
      //              cpl_name,                // Coupling id
      //              1,                       // Postprocessing frequency
      //              CWP_VISU_FORMAT_ENSIGHT, // Postprocessing format
      //              "text");                 // Postprocessing option
    }
  }

  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Create coupling OK\n");
  }



  // Mesh definition
  int          **n_node        = malloc(sizeof(int          *) * n_code);
  double      ***node_coord    = malloc(sizeof(double      **) * n_code);
  int          **n_elt         = malloc(sizeof(int          *) * n_code);
  int         ***elt_node      = malloc(sizeof(int         **) * n_code);
  int          **ijk           = malloc(sizeof(int          *) * n_code);
  PDM_g_num_t ***node_ln_to_gn = malloc(sizeof(PDM_g_num_t **) * n_code);
  PDM_g_num_t ***elt_ln_to_gn  = malloc(sizeof(PDM_g_num_t **) * n_code);

  int *elt_node_idx = NULL;

  for (int i_code = 0; i_code < n_code; i_code++) {
    PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &intra_comm[i_code]);

    double xmin = -0.5;
    double ymin = -0.5;
    double zmin = -0.5;

    double length = 1.;

    if (code_id[i_code] == 2) {
      xmin += 5e-2;
      ymin += 5e-2;
      if (elt_dim == 3) {
        zmin += 5e-2;
      }
      length *= 0.9;
    }

    _gen_mesh(mesh_comm,
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
              &n_node       [i_code],
              &node_coord   [i_code],
              &node_ln_to_gn[i_code],
              &n_elt        [i_code],
              &elt_node     [i_code],
              &elt_ln_to_gn [i_code]);

    int elt_node_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
    if (0) {
      _ijk_grid(elt_type,
                order,
                &ijk[i_code]);
    }
    else {
      int *user_to_ijk = PDM_ho_ordering_user_to_ijk_get("PDM_HO_ORDERING_VTK",
                                                         elt_type,
                                                         order);
      int size = elt_node_n * PDM_Mesh_nodal_elt_dim_get(elt_type);
      ijk[i_code] = malloc(sizeof(int) * size);
      memcpy(ijk[i_code], user_to_ijk, sizeof(int) * size);
    }

    CWP_Block_t     block_type;
    cwipi_element_t old_block_type;
    switch (elt_type) {
    case PDM_MESH_NODAL_POINT:
      block_type     = CWP_BLOCK_NODE;
      old_block_type = CWIPI_NODE;
      break;
    case PDM_MESH_NODAL_BAR2:
      block_type     = CWP_BLOCK_EDGE2;
      old_block_type = CWIPI_EDGE2;
      break;
    case PDM_MESH_NODAL_TRIA3:
      block_type     = CWP_BLOCK_FACE_TRIA3;
      old_block_type = CWIPI_FACE_TRIA3;
      break;
    case PDM_MESH_NODAL_QUAD4:
      block_type     = CWP_BLOCK_FACE_QUAD4;
      old_block_type = CWIPI_FACE_QUAD4;
      break;
    case PDM_MESH_NODAL_POLY_2D:
      block_type     = CWP_BLOCK_FACE_POLY;
      old_block_type = CWIPI_FACE_POLY;
      break;
    case PDM_MESH_NODAL_TETRA4:
      block_type     = CWP_BLOCK_CELL_TETRA4;
      old_block_type = CWIPI_CELL_TETRA4;
      break;
    case PDM_MESH_NODAL_HEXA8:
      block_type     = CWP_BLOCK_CELL_HEXA8;
      old_block_type = CWIPI_CELL_HEXA8;
      break;
    case PDM_MESH_NODAL_PRISM6:
      block_type     = CWP_BLOCK_CELL_PRISM6;
      old_block_type = CWIPI_CELL_PRISM6;
      break;
    case PDM_MESH_NODAL_PYRAMID5:
      block_type     = CWP_BLOCK_CELL_PYRAM5;
      old_block_type = CWIPI_CELL_PYRAM5;
      break;
    case PDM_MESH_NODAL_POLY_3D:
      block_type     = CWP_BLOCK_CELL_POLY;
      old_block_type = CWIPI_NODE;
      break;
    case PDM_MESH_NODAL_BARHO:
      block_type     = CWP_BLOCK_EDGEHO;
      old_block_type = CWIPI_EDGEHO;
      break;
    case PDM_MESH_NODAL_TRIAHO:
      block_type     = CWP_BLOCK_FACE_TRIAHO;
      old_block_type = CWIPI_FACE_TRIAHO;
      break;
    case PDM_MESH_NODAL_QUADHO:
      block_type     = CWP_BLOCK_FACE_QUADHO;
      old_block_type = CWIPI_FACE_QUADHO;
      break;
    case PDM_MESH_NODAL_TETRAHO:
      block_type     = CWP_BLOCK_CELL_TETRAHO;
      old_block_type = CWIPI_CELL_TETRAHO;
      break;
    case PDM_MESH_NODAL_HEXAHO:
      block_type     = CWP_BLOCK_CELL_HEXAHO;
      old_block_type = CWIPI_CELL_HEXAHO;
      break;
    case PDM_MESH_NODAL_PRISMHO:
      block_type     = CWP_BLOCK_CELL_PRISMHO;
      old_block_type = CWIPI_CELL_PRISMHO;
      break;
    case PDM_MESH_NODAL_PYRAMIDHO:
      block_type     = CWP_BLOCK_CELL_PYRAMHO;
      old_block_type = CWIPI_CELL_PYRAMHO;
      break;
    default:
      PDM_error (__FILE__, __LINE__, 0, "invalid elt type %d\n", (int) elt_type);
    }


    if (version == CWP_VERSION_OLD) {

      elt_node_idx = PDM_array_new_idx_from_const_stride_int(elt_node_n, n_elt[i_code][0]);

      cwipi_ho_define_mesh(cpl_name,
                           n_node[i_code][0],
                           n_elt [i_code][0],
                           order,
                           node_coord[i_code][0],
                           elt_node_idx,
                           elt_node[i_code][0]);

      cwipi_ho_ordering_from_IJK_set(cpl_name,
                                     old_block_type,
                                     elt_node_n,
                                     ijk[0]);
    }
    else {
      int block_id = CWP_Mesh_interf_block_add(code_name[i_code],
                                               cpl_name,
                                               block_type);

      for (int i = 0; i < n_part; i++) {
        CWP_Mesh_interf_vtx_set(code_name[i_code],
                                cpl_name,
                                i,
                                n_node       [i_code][i],
                                node_coord   [i_code][i],
                                node_ln_to_gn[i_code][i]);

        CWP_Mesh_interf_block_ho_set(code_name[i_code],
                                     cpl_name,
                                     i,
                                     block_id,
                                     n_elt[i_code][i],
                                     order,
                                     elt_node    [i_code][i],
                                     elt_ln_to_gn[i_code][i]);
      }

      CWP_Mesh_interf_ho_ordering_from_IJK_set(code_name[i_code],
                                               cpl_name,
                                               block_type,
                                               order,
                                               elt_node_n,
                                               ijk[i_code]);

      CWP_Mesh_interf_finalize(code_name[i_code], cpl_name);
    }
  }

  MPI_Barrier(comm);

  if (i_rank == 0) {
    printf("Set mesh OK\n");
  }

  // Create and set fields
  CWP_Status_t visu_status = CWP_STATUS_OFF;
  const char *field_name1 = "field1";

  double ***send_val = malloc(sizeof(double **) * n_code);
  double ***recv_val = malloc(sizeof(double **) * n_code);
  for (int i_code = 0; i_code < n_code; i_code++) {
    send_val[i_code] = malloc(sizeof(double *) * n_part);
    recv_val[i_code] = malloc(sizeof(double *) * n_part);
    for (int i = 0; i < n_part; i++) {
      send_val[i_code][i] = malloc(sizeof(double) * n_node[i_code][i]);
      recv_val[i_code][i] = malloc(sizeof(double) * n_node[i_code][i]);
    }

    if (code_id[i_code] == 1) {
      for (int ipart = 0; ipart < n_part; ipart++) {
        for (int i = 0 ; i < n_node[i_code][ipart]; i++) {
          // send_val[i_code][ipart][i] = node_coord[i_code][ipart][3*i];
          send_val[i_code][ipart][i] = _eval_field(node_coord[i_code][ipart][3*i+0],
                                                   node_coord[i_code][ipart][3*i+1],
                                                   node_coord[i_code][ipart][3*i+2],
                                                   1);
        }
      }
    }
  }

  if (version == CWP_VERSION_NEW) {
    for (int i_code = 0; i_code < n_code; i_code++) {
      if (code_id[i_code] == 1) {
        CWP_Field_create(code_name[i_code],
                         cpl_name,
                         field_name1,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_SEND,
                         visu_status);

        for (int i = 0; i < n_part; i++) {
          CWP_Field_data_set(code_name[i_code],
                             cpl_name,
                             field_name1,
                             i,
                             CWP_FIELD_MAP_SOURCE,
                             send_val[i_code][i]);
        }
      }

      if (code_id[i_code] == 2) {
        CWP_Field_create(code_name[i_code],
                         cpl_name,
                         field_name1,
                         CWP_DOUBLE,
                         CWP_FIELD_STORAGE_INTERLACED,
                         1,
                         CWP_DOF_LOCATION_NODE,
                         CWP_FIELD_EXCH_RECV,
                         visu_status);

        for (int i = 0; i < n_part; i++) {
          CWP_Field_data_set(code_name[i_code],
                             cpl_name,
                             field_name1,
                             i,
                             CWP_FIELD_MAP_TARGET,
                             recv_val[i_code][i]);
        }
      }
    }
  }

  MPI_Barrier(comm);

  if (i_rank == 0) {
    printf("Create fields OK\n");
  }


  PDM_timer_t *timer = PDM_timer_create();
  double t_start, t_end;
  MPI_Barrier(comm);

  PDM_timer_init(timer);

  t_start = PDM_timer_elapsed(timer);
  PDM_timer_resume(timer);

  if (version == CWP_VERSION_OLD) {
    cwipi_locate(cpl_name);
  }
  else {
    char str_tolerance[99];
    sprintf(str_tolerance, "%f", tolerance);
    for (int i_code = 0; i_code < n_code; i_code++) {
      CWP_Spatial_interp_property_set(code_name[i_code], cpl_name, "tolerance", CWP_DOUBLE, str_tolerance);
    }

    for (int i_code = 0; i_code < n_code; i_code++) {
      CWP_Spatial_interp_weights_compute(code_name[i_code], cpl_name);
    }
  }

  MPI_Barrier(comm);

  PDM_timer_hang_on(timer);
  t_end = PDM_timer_elapsed(timer);

  PDM_timer_free(timer);

  double dt = t_end - t_start;
  double dt_min;
  double dt_max;
  MPI_Reduce(&dt, &dt_min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
  MPI_Reduce(&dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);


  MPI_Barrier(comm);
  if (i_rank == 0) {
    printf("Interpolation weights computation OK\n");
    printf("dt_min = %12.5es\n", dt_min);
    printf("dt_max = %12.5es\n", dt_max);
  }

  MPI_Barrier(comm);

  if (version == CWP_VERSION_OLD) {
    int request;
    if (code_id[0] == 1) {
      cwipi_issend(cpl_name, "ech", 0, 1, 1, 0.1, field_name1, send_val[0][0], &request);
    }
    else {
      cwipi_irecv (cpl_name, "ech", 0, 1, 1, 0.1, field_name1, recv_val[0][0], &request);
    }

    if (code_id[0] == 1) {
      cwipi_wait_issend(cpl_name, request);
    }
    else {
      cwipi_wait_irecv (cpl_name, request);
    }
  }
  else {
    for (int i_code = 0; i_code < n_code; i_code++) {
      if (code_id[i_code] == 1) {
        CWP_Field_issend(code_name[i_code], cpl_name, field_name1);
      }
      if (code_id[i_code] == 2) {
        CWP_Field_irecv (code_name[i_code], cpl_name, field_name1);
      }


      if (code_id[i_code] == 1) {
        CWP_Field_wait_issend(code_name[i_code], cpl_name, field_name1);
      }
      if (code_id[i_code] == 2) {
        CWP_Field_wait_irecv (code_name[i_code], cpl_name, field_name1);
      }
    }
  }

  if (i_rank == 0) {
    printf("Exchange fields OK\n");
  }

  //  Check
  PDM_g_num_t gn_wrong = 0;


  if (version == CWP_VERSION_OLD ||
      (spatial_interp == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT ||
       spatial_interp == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE         ||
       spatial_interp == CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE)) {
    double max_err = 0.;
    PDM_g_num_t n_wrong = 0;

    for (int i_code = 0; i_code < n_code; i_code++) {
      if (code_id[i_code] == 2) {
        for (int ipart = 0; ipart < n_part; ipart++) {
          int n_computed_tgt = 0;
          const int *computed_tgt = NULL;

          if (version == CWP_VERSION_OLD) {
            n_computed_tgt = cwipi_get_n_located_points(cpl_name);
            computed_tgt   = cwipi_get_located_points  (cpl_name);
          }
          else {
            n_computed_tgt = CWP_N_computed_tgts_get(code_name[i_code],
                                                     cpl_name,
                                                     field_name1,
                                                     ipart);
            computed_tgt = CWP_Computed_tgts_get(code_name[i_code],
                                                 cpl_name,
                                                 field_name1,
                                                 ipart);
          }

          // log_trace("part %d, n_computed_tgt = %d / %d\n", ipart, n_computed_tgt, n_node[i_code][ipart]);

          for (int i = 0; i < n_computed_tgt; i++) {
            int j = computed_tgt[i] - 1;

            double val = _eval_field(node_coord[i_code][ipart][3*i+0],
                                     node_coord[i_code][ipart][3*i+1],
                                     node_coord[i_code][ipart][3*i+2],
                                     1);

            double err = ABS(recv_val[i_code][ipart][j] - val);

            max_err = MAX(max_err, err);

            n_wrong += (err > 1e-6);
          }
        }
      }
    }

    PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) &comm);
    PDM_MPI_Allreduce(&n_wrong, &gn_wrong, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, pdm_comm);

    double gmax_err;
    MPI_Allreduce(&max_err, &gmax_err, 1, MPI_DOUBLE, MPI_MAX, comm);

    if (i_rank == 0) {
      printf("Max error = %g ("PDM_FMT_G_NUM" wrong points)\n", gmax_err, gn_wrong);
    }
  }


  if (0) {
    // Visu VTK
    for (int i_code = 0; i_code < n_code; i_code++) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        char filename[999];
        sprintf(filename, "c_new_vs_old_HO_%s_%d_%d.vtk", code_name[i_code], i_part, i_rank);

        double *field_value;
        if (code_id[i_code] == 1) {
          field_value = send_val[i_code][i_part];
        }
        else {
          field_value = recv_val[i_code][i_part];
        }

        PDM_vtk_write_std_elements_ho_with_vtx_field(filename,
                                                     order,
                                                     n_node       [i_code][i_part],
                                                     node_coord   [i_code][i_part],
                                                     node_ln_to_gn[i_code][i_part],
                                                     elt_type,
                                                     n_elt       [i_code][i_part],
                                                     elt_node    [i_code][i_part],
                                                     elt_ln_to_gn[i_code][i_part],
                                                     0,
                                                     NULL,
                                                     NULL,
                                                     1,
                                                     &field_name1,
                                   (const double **) &field_value);
      }
    }
  }



  for (int i_code = 0; i_code < n_code; i_code++) {
    if (version == CWP_VERSION_OLD) {
      cwipi_delete_coupling(cpl_name);

      cwipi_finalize();
    }
    else {
      CWP_Mesh_interf_del(code_name[i_code], cpl_name);

      CWP_Cpl_del(code_name[i_code], cpl_name);

      CWP_Finalize();
    }

    for (int i = 0; i < n_part; i++) {
      free(node_coord   [i_code][i]);
      free(node_ln_to_gn[i_code][i]);
      free(elt_node     [i_code][i]);
      free(elt_ln_to_gn [i_code][i]);
      free(send_val     [i_code][i]);
      free(recv_val     [i_code][i]);
    }
    free(n_node       [i_code]);
    free(node_coord   [i_code]);
    free(node_ln_to_gn[i_code]);
    free(n_elt        [i_code]);
    free(elt_node     [i_code]);
    free(elt_ln_to_gn [i_code]);
    free(send_val     [i_code]);
    free(recv_val     [i_code]);
    free(ijk          [i_code]);
  }
  free(n_node       );
  free(node_coord   );
  free(node_ln_to_gn);
  free(n_elt        );
  free(elt_node     );
  free(elt_ln_to_gn );
  free(send_val     );
  free(recv_val     );
  free(ijk          );

  if (version == CWP_VERSION_OLD) {
    free(elt_node_idx);
  }

  free(coupled_code_name);
  free(code_name);
  free(intra_comm);

  // Finalize MPI
  MPI_Finalize();

  return (gn_wrong > 0);
}
