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




static double
_eval_field
(
 const double  x,
 const double  y,
 const double  z
 )
{
  const double x0 = 0.;
  const double y0 = 0.;
  const double z0 = 0.;
  const double r0 = 0.48;

  const double x1 = -0.18;
  const double y1 = 0.18;
  const double z1 = 0.32;
  const double r1 = 0.11;

  const double x2 = -x1;
  const double y2 = y1;
  const double z2 = z1;
  const double r2 = r1;

  const double fm1 = -0.2*z + y + 0.18;
  const double fm2 = -(y + z) + 0.13;
  double fm;
  if (fm1 > fm2) {
    fm = -fm1;
  } else {
    fm = -fm2;
  }

  double F[4] = {
    (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) - r0*r0,
    -((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1) - r1*r1),
    -((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2) - r2*r2),
    fm
  };

  double f = -HUGE_VAL;
  for (int i = 0; i < 4; i++) {
    f = PDM_MAX(f, F[i]);
  }

  return f;
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
int                        **edge_vtx,
double                     **vtx_coord,
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



  PDM_g_num_t *distrib_vtx = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];

  *n_vtx = dn_vtx;
  *vtx_coord = malloc(sizeof(double) * dn_vtx * 3);
  memcpy(*vtx_coord, dvtx_coord, sizeof(double) * dn_vtx * 3);



  PDM_dmesh_nodal_elmts_t *dmne = dmn->volumic;
  int id_section = dmne->sections_id[0];

  PDM_g_num_t *dcell_vtx = PDM_DMesh_nodal_elmts_section_std_get(dmne, id_section);

  int cell_vtx_n = PDM_Mesh_nodal_n_vertices_element(elt_type, 1);
  *cell_vtx = malloc(sizeof(int) * dn_cell * cell_vtx_n);
  for (int i = 0; i < dn_cell * cell_vtx_n; i++) {
    (*cell_vtx)[i] = (int) dcell_vtx[i];
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
   * ./pdm_t_tag_cells_sol -n <n_vtx_seg> -t <elt_type>
   *
   * elt_type :
   *  0 -> tetrahedra
   *  1 -> pyramids
   *  2 -> prisms
   *  3 -> hexahedra
   */

  PDM_g_num_t n_vtx_seg = 30; // Number of vtx on each side of the cube mesh
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

  int     n_cell        = 0;
  int     n_face        = 0;
  int     n_edge        = 0;
  int     n_vtx         = 0;
  int    *cell_face_idx = NULL;
  int    *cell_face     = NULL;
  int    *face_edge_idx = NULL;
  int    *face_edge     = NULL;
  int    *edge_vtx      = NULL;
  double *vtx_coord     = NULL;
  int    *cell_vtx      = NULL;
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
                 &edge_vtx,
                 &vtx_coord,
                 &cell_vtx);


  /* Compute field values at vertices */
  int *vtx_field = malloc(sizeof(int) * n_vtx);
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    double val = _eval_field(vtx_coord[3*ivtx  ],
                             vtx_coord[3*ivtx+1],
                             vtx_coord[3*ivtx+2]);
    vtx_field[ivtx] = PDM_SIGN(val);
  }

  /*******************************************************************
   *
   * /!\ 'cell_vtx' may only be used for visualization purposes /!\
   *
   *******************************************************************/

  /*
   *  1) Tag cells
   *
   *    cell_tag = 1 if some of the cell's vtx have positive values (●)
   *                 and some have negative values (○)
   *               0 else
   *
   *    ● ‒ ● ‒ ● ‒ ● ‒ ○ ‒ ○ ‒ ○
   *    | 0 | 0 | 0 | 1 | 0 | 0 |
   *    ● ‒ ● ‒ ● ‒ ● ‒ ○ ‒ ○ ‒ ○
   *    | 0 | 0 | 1 | 1 | 0 | 0 |
   *    ● ‒ ● ‒ ● ‒ ○ ‒ ○ ‒ ○ ‒ ○
   *    | 0 | 1 | 1 | 0 | 0 | 0 |
   *    ● ‒ ● ‒ ○ ‒ ○ ‒ ○ ‒ ○ ‒ ○
   *    | 0 | 1 | 0 | 0 | 0 | 0 |
   *    ● ‒ ● ‒ ○ ‒ ○ ‒ ○ ‒ ○ ‒ ○
   *
   */

  int *cell_tag = malloc(sizeof(int) * n_cell);


  for (int icell = 0; icell < n_cell; icell++) {

    cell_tag[icell] = 0;

    for (int idx_face = cell_face_idx[icell]; idx_face < cell_face_idx[icell+1]; idx_face++) {
      int iface = PDM_ABS(cell_face[idx_face]) - 1;

      for (int idx_edge = face_edge_idx[iface]; idx_edge < face_edge_idx[iface+1]; idx_edge++) {
        int iedge = PDM_ABS(face_edge[idx_edge]) - 1;

        int ivtx1 = edge_vtx[2*iedge  ] - 1;
        int ivtx2 = edge_vtx[2*iedge+1] - 1;

        if (PDM_SIGN(vtx_field[ivtx1]) != PDM_SIGN(vtx_field[ivtx2])) {
          cell_tag[icell] = 1;
          break; // break out of loop on current face's edges
        }
      } // end of loop on current face's edges

      if (cell_tag[icell] == 1) {
        break; // break out of loop on current cell's faces
      }
    } // end of loop on current cell's faces
  } // end of loop on cells



  /*
   *  Visualize cell tags
   */
  const char* field_name[] = {"cell_tag", 0 };
  if(0 == 1) {
    PDM_vtk_write_std_elements("cell_tag.vtk",
                               n_vtx,
                               vtx_coord,
                               NULL,
                               elt_type,
                               n_cell,
                               cell_vtx,
                               NULL,
                               1,
                               field_name,
                (const int **) &cell_tag);
  }




  /*
   *  2) Tag all cells
   *
   *    cells sharing a face with a cell tagged 'n' get the tag 'n+1'
   *
   *    ● ‒ ● ‒ ● ‒ ● ‒ ○ ‒ ○ ‒ ○
   *    | 4 | 3 | 2 | 1 | 2 | 3 |
   *    ● ‒ ● ‒ ● ‒ ● ‒ ○ ‒ ○ ‒ ○
   *    | 3 | 2 | 1 | 1 | 2 | 3 |
   *    ● ‒ ● ‒ ● ‒ ○ ‒ ○ ‒ ○ ‒ ○
   *    | 2 | 1 | 1 | 2 | 3 | 4 |
   *    ● ‒ ● ‒ ○ ‒ ○ ‒ ○ ‒ ○ ‒ ○
   *    | 2 | 1 | 2 | 3 | 4 | 5 |
   *    ● ‒ ● ‒ ○ ‒ ○ ‒ ○ ‒ ○ ‒ ○
   *
   */

  int *face_cell_idx = NULL;
  int *face_cell     = NULL;
  PDM_connectivity_transpose(n_cell,
                             n_face,
                             cell_face_idx,
                             cell_face,
                             &face_cell_idx,
                             &face_cell);

  int *newly_tagged_cells = malloc(sizeof(int) * n_cell);
  int n_newly_tagged_cells = 0;
  for (int icell = 0; icell < n_cell; icell++) {
    if (cell_tag[icell] != 0) {
      newly_tagged_cells[n_newly_tagged_cells++] = icell;
    }
  }

  int *tmp_cells = malloc(sizeof(int) * n_cell);


  int i_step = 1;
  while (n_newly_tagged_cells > 0) {

    printf("step = %d, n = %d\n", i_step, n_newly_tagged_cells);

    int tmp_n = 0;

    for (int i = 0; i < n_newly_tagged_cells; i++) {

      int icell = newly_tagged_cells[i];

      for (int idx_face = cell_face_idx[icell]; idx_face < cell_face_idx[icell+1]; idx_face++) {
        int iface = PDM_ABS(cell_face[idx_face]) - 1;

        for (int idx_cell = face_cell_idx[iface]; idx_cell < face_cell_idx[iface+1]; idx_cell++) {
          int jcell = PDM_ABS(face_cell[idx_cell]) - 1;

          if (cell_tag[jcell] == 0) {
            tmp_cells[tmp_n++] = jcell;
            cell_tag[jcell] = i_step + 1;
          }
        } // end of loop on current face's cells

      } // end of loop on current cell's faces
    } // end of loop on cells

    memcpy(newly_tagged_cells, tmp_cells, sizeof(int) * tmp_n);
    n_newly_tagged_cells = tmp_n;

    i_step++;
  }
  free(tmp_cells);
  free(newly_tagged_cells);
  free(face_cell_idx);
  free(face_cell);


  /*
   *  Visualize 'expanded' cell tags
   */
  if(1 == 0) {
    PDM_vtk_write_std_elements("expanded_cell_tag.vtk",
                               n_vtx,
                               vtx_coord,
                               NULL,
                               elt_type,
                               n_cell,
                               cell_vtx,
                               NULL,
                               1,
                               field_name,
                (const int **) &cell_tag);
  }


  /* Free memory */
  free(cell_face_idx);
  free(cell_face    );
  free(face_edge_idx);
  free(face_edge    );
  free(edge_vtx     );
  free(vtx_coord    );
  free(cell_vtx     );
  free(vtx_field    );
  free(cell_tag     );

  PDM_MPI_Finalize();
  return 0;
}


