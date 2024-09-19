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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_dmesh.h"
#include "pdm_part_connectivity_transform.h"

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
     "  -n      <n_vtx_seg> Number of vertices on each side of the cube mesh.\n\n"
     "  -t      <elt_type>  Type of cells.\n\n"
     "  -h                  This message.\n\n");

  exit(exit_code);
}

static void
_read_args(int           argc,
           char        **argv,
           PDM_g_num_t  *n_vtx_seg,
           int          *elt_type)
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
        *n_vtx_seg = (PDM_g_num_t) atol(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = atoi(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_compute_face_vtx
(
 const int   n_face,
 int        *face_edge_idx,
 int        *face_edge,
 int        *edge_vtx,
 int       **face_vtx
 )
{
  int dbg = 0;
  *face_vtx = (int *) malloc(sizeof(int) * face_edge_idx[n_face]);

  for (int i = 0; i < n_face; i++) {

    if (dbg) {
      log_trace("\nFace %d\n", i);
      for (int idx_edge = face_edge_idx[i]; idx_edge < face_edge_idx[i+1]; idx_edge++) {
        int iedge = PDM_ABS(face_edge[idx_edge]) - 1;
        log_trace("  edge %d: %d %d\n",
                  face_edge[idx_edge],
                  edge_vtx[2*iedge], edge_vtx[2*iedge+1]);
      }
    }
    int *fv = *face_vtx + face_edge_idx[i];

    int cur_vtx, next_vtx;
    int cur_edge = face_edge[face_edge_idx[i]];
    if (cur_edge < 0) {
      cur_edge = -cur_edge - 1;
      cur_vtx  = edge_vtx[2*cur_edge+1];
      next_vtx = edge_vtx[2*cur_edge  ];
    } else {
      cur_edge = cur_edge - 1;
      cur_vtx  = edge_vtx[2*cur_edge  ];
      next_vtx = edge_vtx[2*cur_edge+1];
    }

    for (int ivtx = 0; ivtx < face_edge_idx[i+1] - face_edge_idx[i]; ivtx++) {
      fv[ivtx] = cur_vtx;

      for (int iedg = face_edge_idx[i]; iedg < face_edge_idx[i+1]; iedg++) {
        cur_edge = face_edge[iedg];
        int vtx1, vtx2;
        if (cur_edge < 0) {
          cur_edge = -cur_edge - 1;
          vtx1 = edge_vtx[2*cur_edge+1];
          vtx2 = edge_vtx[2*cur_edge  ];
        } else {
          cur_edge = cur_edge - 1;
          vtx1 = edge_vtx[2*cur_edge  ];
          vtx2 = edge_vtx[2*cur_edge+1];
        }

        if (vtx1 == next_vtx) {
          cur_vtx  = next_vtx;
          next_vtx = vtx2;
          break;
        }
      }
    }

    if (dbg) {
      log_trace("  face_vtx = ");
      for (int ivtx = 0; ivtx < face_edge_idx[i+1] - face_edge_idx[i]; ivtx++) {
        log_trace("%d ", fv[ivtx]);
      }
      log_trace("\n");
    }

  }
}


static void
_generate_mesh
(
 const PDM_g_num_t            n_vtx_seg,
 const PDM_Mesh_nodal_elt_t   elt_type,
 int                         *n_cell,
 int                         *n_face,
 int                         *n_edge,
 int                         *n_vtx,
 int                        **cell_face_idx,
 int                        **cell_face,
 int                        **face_edge_idx,
 int                        **face_edge,
 int                        **edge_vtx_idx,
 int                        **edge_vtx,
 double                     **vtx_coord,
 int                         *n_face_group,
 int                        **group_face_idx,
 int                        **group_face,
 int                         *n_edge_group,
 int                        **group_edge_idx,
 int                        **group_edge,
 int                        **cell_vtx
)
{
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         1,
                                                         -0.5,
                                                         -0.5,
                                                         -0.5,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_dmesh_nodal_to_dmesh_t *dmntodm = PDM_dmesh_nodal_to_dmesh_create(1,
                                                                        comm,
                                                                        PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  PDM_dmesh_t *dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);


  int          dn_cell = 0;
  int          dn_face = 0;
  int          dn_edge = 0;
  int          dn_vtx  = 0;
  double      *dvtx_coord     = NULL;
  int         *dcell_face_idx = NULL;
  PDM_g_num_t *dcell_face     = NULL;
  int         *dface_edge_idx = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  int         *dedge_vtx_idx  = NULL;
  PDM_g_num_t *dedge_vtx      = NULL;

  dn_cell = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                       &dcell_face,
                                       &dcell_face_idx,
                                       PDM_OWNERSHIP_KEEP);

  dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                       &dface_edge,
                                       &dface_edge_idx,
                                       PDM_OWNERSHIP_KEEP);

  dn_edge = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                       &dedge_vtx,
                                       &dedge_vtx_idx,
                                       PDM_OWNERSHIP_KEEP);


  *n_cell = dn_cell;
  *cell_face_idx = malloc(sizeof(int) * (dn_cell + 1));
  memcpy(*cell_face_idx, dcell_face_idx, sizeof(int) * (dn_cell + 1));
  *cell_face = malloc(sizeof(int) * dcell_face_idx[dn_cell]);
  for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
    (*cell_face)[i] = (int) dcell_face[i];
  }


  *n_face = dn_face;
  *face_edge_idx = malloc(sizeof(int) * (dn_face + 1));
  memcpy(*face_edge_idx, dface_edge_idx, sizeof(int) * (dn_face + 1));
  *face_edge = malloc(sizeof(int) * dface_edge_idx[dn_face]);
  for (int i = 0; i < dface_edge_idx[dn_face]; i++) {
    (*face_edge)[i] = (int) dface_edge[i];
  }


  *n_edge = dn_edge;
  *edge_vtx = malloc(sizeof(int) * 2*dn_edge);
  for (int i = 0; i < 2*dn_edge; i++) {
    (*edge_vtx)[i] = (int) dedge_vtx[i];
  }
  *edge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);



  PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];

  *n_vtx = dn_vtx;
  *vtx_coord = malloc(sizeof(double) * dn_vtx * 3);
  memcpy(*vtx_coord, dvtx_coord, sizeof(double) * dn_vtx * 3);


  /* cell -> vtx */
  PDM_dmesh_nodal_elmts_t *dmne = dmn->volumic;
  int id_section = dmne->sections_id[0];

  PDM_g_num_t *dcell_vtx = PDM_DMesh_nodal_elmts_section_std_get(dmne, id_section);

  int cell_vtx_n = PDM_Mesh_nodal_n_vertices_element(elt_type, 1);
  *cell_vtx = malloc(sizeof(int) * dn_cell * cell_vtx_n);
  for (int i = 0; i < dn_cell * cell_vtx_n; i++) {
    (*cell_vtx)[i] = (int) dcell_vtx[i];
  }



  /* Groups */
  int         *dgroup_face_idx = NULL;
  PDM_g_num_t *dgroup_face     = NULL;
  *n_face_group = PDM_dmesh_bound_get(dmesh,
                                     PDM_BOUND_TYPE_FACE,
                                     &dgroup_face,
                                     &dgroup_face_idx,
                                     PDM_OWNERSHIP_KEEP);

  int         *dgroup_edge_idx = NULL;
  PDM_g_num_t *dgroup_edge     = NULL;
  *n_edge_group = PDM_dmesh_bound_get(dmesh,
                                     PDM_BOUND_TYPE_EDGE,
                                     &dgroup_edge,
                                     &dgroup_edge_idx,
                                     PDM_OWNERSHIP_KEEP);


  *group_face_idx = malloc(sizeof(int) * (*n_face_group + 1));
  memcpy(*group_face_idx, dgroup_face_idx, sizeof(int) * (*n_face_group + 1));
  *group_face = malloc(sizeof(int) * dgroup_face_idx[*n_face_group]);
  for (int i = 0; i < dgroup_face_idx[*n_face_group]; i++) {
    (*group_face)[i] = (int) dgroup_face[i];
  }


  *group_edge_idx = malloc(sizeof(int) * (*n_edge_group + 1));
  memcpy(*group_edge_idx, dgroup_edge_idx, sizeof(int) * (*n_edge_group + 1));
  *group_edge = malloc(sizeof(int) * dgroup_edge_idx[*n_edge_group]);
  for (int i = 0; i < dgroup_edge_idx[*n_edge_group]; i++) {
    (*group_edge)[i] = (int) dgroup_edge[i];
  }


  PDM_dcube_nodal_gen_free(dcube);
  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
}




/**
 *
 * \brief  Main
 *
 */
int main(int argc, char *argv[])
{

  /*
   * ./pdm_t_groups_sol -n <n_vtx_seg> -t <elt_type>
   *
   * elt_type :
   *  0 -> tetrahedra
   *  1 -> pyramids
   *  2 -> prisms
   *  3 -> hexahedra
   */

  PDM_g_num_t n_vtx_seg = 10; // Number of vtx on each side of the cube mesh
  int         _elt_type = 0;  // Type of cells
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &_elt_type);


  PDM_MPI_Init (&argc, &argv);

  PDM_Mesh_nodal_elt_t elt_type;
  switch (_elt_type  ) {
    case 0:
    elt_type = PDM_MESH_NODAL_TETRA4;
    break;

    case 1:
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    break;

    case 2:
    elt_type = PDM_MESH_NODAL_PRISM6;
    break;

    default:
    elt_type = PDM_MESH_NODAL_HEXA8;
  }

  int     n_cell         = 0;
  int     n_face         = 0;
  int     n_edge         = 0;
  int     n_vtx          = 0;
  int    *cell_face_idx  = NULL;
  int    *cell_face      = NULL;
  int    *face_edge_idx  = NULL;
  int    *face_edge      = NULL;
  int    *edge_vtx_idx   = NULL;
  int    *edge_vtx       = NULL;
  double *vtx_coord      = NULL;
  int     n_face_group   = 0;
  int    *group_face_idx = NULL;
  int    *group_face     = NULL;
  int     n_edge_group   = 0;
  int    *group_edge_idx = NULL;
  int    *group_edge     = NULL;
  int    *cell_vtx       = NULL;
  _generate_mesh(n_vtx_seg,
                 elt_type,
                 &n_cell,
                 &n_face,
                 &n_edge,
                 &n_vtx,
                 &cell_face_idx,
                 &cell_face,
                 &face_edge_idx,
                 &face_edge,
                 &edge_vtx_idx,
                 &edge_vtx,
                 &vtx_coord,
                 &n_face_group,
                 &group_face_idx,
                 &group_face,
                 &n_edge_group,
                 &group_edge_idx,
                 &group_edge,
                 &cell_vtx);


  /*
   *  Goals :
   *    - transfer face groups to edges and vertices
   *    - transfer edge groups to vertices
   *
   */

  int *edge_surface_idx = NULL;
  int *edge_surface     = NULL;
  int *vtx_surface_idx  = NULL;
  int *vtx_surface      = NULL;
  int *vtx_ridge_idx    = NULL;
  int *vtx_ridge        = NULL;


  /* Transpose group_face to get face_group */
  int *face_group_idx = NULL;
  int *face_group     = NULL;
  PDM_connectivity_transpose(n_face_group,
                             n_face,
                             group_face_idx,
                             group_face,
                             &face_group_idx,
                             &face_group);


  /* Transpose face_edge to get edge_face */
  int *edge_face_idx = NULL;
  int *edge_face     = NULL;
  PDM_connectivity_transpose(n_face,
                             n_edge,
                             face_edge_idx,
                             face_edge,
                             &edge_face_idx,
                             &edge_face);


  /* Combine edge_face and face_group to get edge_surface  */
  PDM_combine_connectivity(n_edge,
                           edge_face_idx,
                           edge_face,
                           face_group_idx,
                           face_group,
                           &edge_surface_idx,
                           &edge_surface);


  /* Transpose group_edge to get edge_group */
  int *edge_group_idx = NULL;
  int *edge_group     = NULL;
  PDM_connectivity_transpose(n_edge_group,
                             n_edge,
                             group_edge_idx,
                             group_edge,
                             &edge_group_idx,
                             &edge_group);

  /* Transpose edge_vtx to get vtx_edge */
  int *vtx_edge_idx = NULL;
  int *vtx_edge     = NULL;
  PDM_connectivity_transpose(n_edge,
                             n_vtx,
                             edge_vtx_idx,
                             edge_vtx,
                             &vtx_edge_idx,
                             &vtx_edge);

  /* Combine vtx_edge and edge_group to get vtx_ridge  */
  PDM_combine_connectivity(n_vtx,
                           vtx_edge_idx,
                           vtx_edge,
                           edge_group_idx,
                           edge_group,
                           &vtx_ridge_idx,
                           &vtx_ridge);


  /* Combine vtx_edge and edge_surface to get vtx_surface  */
  PDM_combine_connectivity(n_vtx,
                           vtx_edge_idx,
                           vtx_edge,
                           edge_surface_idx,
                           edge_surface,
                           &vtx_surface_idx,
                           &vtx_surface);

  if(1 == 0) {
    PDM_log_trace_connectivity_int(edge_surface_idx,
                                   edge_surface,
                                   n_edge,
                                   "edge_surface : ");

    PDM_log_trace_connectivity_int(vtx_ridge_idx,
                                   vtx_ridge,
                                   n_vtx,
                                   "vtx_ridge : ");

    PDM_log_trace_connectivity_int(vtx_surface_idx,
                                   vtx_surface,
                                   n_vtx,
                                   "vtx_surface : ");
  }

  /* Visu */
  int *face_vtx = NULL;
  _compute_face_vtx(n_face,
                    face_edge_idx,
                    face_edge,
                    edge_vtx,
                    &face_vtx);

  int *face_group_n = (int *) malloc(sizeof(int) * n_face);
  for (int i = 0; i < n_face; i++) {
    face_group_n[i] = face_group_idx[i+1] - face_group_idx[i];
  }

  int *edge_group_n   = (int *) malloc(sizeof(int) * n_edge);
  int *edge_surface_n = (int *) malloc(sizeof(int) * n_edge);
  for (int i = 0; i < n_edge; i++) {
    edge_group_n[i] = edge_group_idx[i+1] - edge_group_idx[i];
    edge_surface_n[i] = edge_surface_idx[i+1] - edge_surface_idx[i];
  }

  int *vtx_ridge_n   = (int *) malloc(sizeof(int) * n_vtx);
  int *vtx_surface_n = (int *) malloc(sizeof(int) * n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    vtx_ridge_n[i]   = vtx_ridge_idx[i+1]   - vtx_ridge_idx[i];
    vtx_surface_n[i] = vtx_surface_idx[i+1] - vtx_surface_idx[i];
  }

  if(0 == 1) {
    PDM_vtk_write_polydata("face_group.vtk",
                           n_vtx,
                           vtx_coord,
                           NULL,
                           n_face,
                           face_edge_idx,
                           face_vtx,
                           NULL,
                           face_group_n);

    const char* field_name[] = {"n_ridge", "n_surface", 0 };

    const int *edge_field[2] = {edge_group_n, edge_surface_n};
    PDM_vtk_write_std_elements("edge_group.vtk",
                               n_vtx,
                               vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_BAR2,
                               n_edge,
                               edge_vtx,
                               NULL,
                               2,
                               field_name,
                               edge_field);

    const int *vtx_field[2] = {vtx_ridge_n, vtx_surface_n};
    PDM_vtk_write_std_elements("vtx_group.vtk",
                               n_vtx,
                               vtx_coord,
                               NULL,
                               PDM_MESH_NODAL_POINT,
                               n_vtx,
                               NULL,
                               NULL,
                               2,
                               field_name,
                               vtx_field);
  }
  free(face_vtx);
  free(face_group_n);
  free(edge_group_n);
  free(edge_surface_n);
  free(vtx_ridge_n);
  free(vtx_surface_n);

  free(face_group_idx);
  free(face_group);
  free(edge_group_idx);
  free(edge_group);
  free(edge_face_idx);
  free(edge_face);
  free(vtx_edge_idx);
  free(vtx_edge);

  free(edge_surface_idx);
  free(edge_surface    );
  free(vtx_surface_idx );
  free(vtx_surface     );
  free(vtx_ridge_idx   );
  free(vtx_ridge       );


  /* Free memory */
  free(cell_face_idx );
  free(cell_face     );
  free(face_edge_idx );
  free(face_edge     );
  free(edge_vtx_idx  );
  free(edge_vtx      );
  free(vtx_coord     );
  free(group_face_idx);
  free(group_face    );
  free(group_edge_idx);
  free(group_edge    );
  free(cell_vtx      );

  PDM_MPI_Finalize();
  return 0;
}
