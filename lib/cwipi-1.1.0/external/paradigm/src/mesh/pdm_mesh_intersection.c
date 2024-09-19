/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mesh_intersection_priv.h"
#include "pdm_mesh_intersection.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_extract_part.h"
#include "pdm_extract_part_priv.h"
#include "pdm_part_to_part.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_writer.h"
#include "pdm_gnum.h"
#include "pdm_box_priv.h"
#include "pdm_unique.h"
#include "pdm_triangulate.h"
#include "pdm_geom_elem.h"
#include "pdm_line.h"
#include "pdm_plane.h"
#include "pdm_polygon.h"
#include "pdm_predicate.h"
#include "pdm_order.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_part_mesh_nodal_to_pmesh.h"

#include "pdm_mesh_intersection_surf_surf_atomic.h"
#include "pdm_mesh_intersection_vol_vol_atomic.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static
void
_export_vtk_1d
(
 const char               *pattern,
       PDM_extract_part_t *extrp_mesh
)
{
  int i_rank;
  PDM_MPI_Comm_rank(extrp_mesh->comm, &i_rank);


  for(int i_part = 0; i_part < extrp_mesh->n_part_out; ++i_part) {

    PDM_g_num_t *edge_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;
    int n_edge = PDM_extract_part_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_EDGE  , &edge_ln_to_gn, PDM_OWNERSHIP_KEEP);
    int n_vtx  = PDM_extract_part_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_VTX, &vtx_ln_to_gn , PDM_OWNERSHIP_KEEP);

    double *vtx_coord = NULL;
    PDM_extract_part_vtx_coord_get(extrp_mesh, i_part, &vtx_coord, PDM_OWNERSHIP_KEEP);

    int  *edge_vtx      = NULL;
    int  *edge_vtx_idx  = NULL;
    PDM_extract_part_connectivity_get(extrp_mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , &edge_vtx , &edge_vtx_idx , PDM_OWNERSHIP_KEEP);

    char filename[999];
    sprintf(filename, "%s_%i_%i.vtk", pattern, i_part, i_rank);
    PDM_vtk_write_std_elements(filename,
                               n_vtx,
                               vtx_coord,
                               vtx_ln_to_gn,
                               PDM_MESH_NODAL_BAR2,
                               n_edge,
                               edge_vtx,
                               edge_ln_to_gn,
                               0,
                               NULL,
                               NULL);
  }
}

static
void
_export_vtk_2d
(
 const char               *pattern,
       PDM_extract_part_t *extrp_mesh
)
{
  int i_rank;
  PDM_MPI_Comm_rank(extrp_mesh->comm, &i_rank);

  for(int i_part = 0; i_part < extrp_mesh->n_part_out; ++i_part) {

    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;
    int n_face = extrp_mesh->n_target[i_part];
    face_ln_to_gn = extrp_mesh->target_gnum[i_part];
    // int n_face = PDM_extract_part_parent_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_FACE  , &face_ln_to_gn, PDM_OWNERSHIP_KEEP);
    int n_vtx  = PDM_extract_part_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_VTX, &vtx_ln_to_gn , PDM_OWNERSHIP_KEEP);

    double *vtx_coord = NULL;
    PDM_extract_part_vtx_coord_get(extrp_mesh, i_part, &vtx_coord, PDM_OWNERSHIP_KEEP);

    int  *face_edge     = NULL;
    int  *face_edge_idx = NULL;
    int  *edge_vtx      = NULL;
    int  *edge_vtx_idx  = NULL;
    PDM_extract_part_connectivity_get(extrp_mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE, &face_edge, &face_edge_idx, PDM_OWNERSHIP_KEEP);
    PDM_extract_part_connectivity_get(extrp_mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , &edge_vtx , &edge_vtx_idx , PDM_OWNERSHIP_KEEP);

    int *face_vtx = NULL;
    PDM_compute_face_vtx_from_face_and_edge(n_face, face_edge_idx, face_edge, edge_vtx, &face_vtx);

    char filename[999];
    sprintf(filename, "%s_%i_%i.vtk", pattern, i_part, i_rank);
    PDM_vtk_write_polydata(filename,
                           n_vtx,
                           vtx_coord,
                           vtx_ln_to_gn,
                           n_face,
                           face_edge_idx,
                           face_vtx,
                           face_ln_to_gn,
                           NULL);


    free(face_vtx);
  }
}

static void
_export_vtk_3d
(
 const char               *name_chr,
       PDM_extract_part_t *extrp
)
{
  int i_rank;

  PDM_MPI_Comm_rank(extrp->comm, &i_rank);

  int          *pn_extract_cell        = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
  int          *pn_extract_face        = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
  int          *pn_extract_vtx         = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
  int         **pextract_cell_face     = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  int         **pextract_cell_face_idx = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  int         **pextract_face_vtx      = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  int         **pextract_face_vtx_idx  = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  double      **pextract_vtx           = (double      **) malloc(extrp->n_part_out * sizeof(double      *));
  PDM_g_num_t **pextract_cell_ln_to_gn = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pextract_face_ln_to_gn = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pextract_vtx_ln_to_gn  = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));


  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {

    pn_extract_cell[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_CELL);

    pn_extract_face[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE);

    pn_extract_vtx[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                           i_part,
                                                           PDM_MESH_ENTITY_VTX);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                      &pextract_cell_face[i_part],
                                      &pextract_cell_face_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                      &pextract_face_vtx[i_part],
                                      &pextract_face_vtx_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    PDM_ownership_t owner = PDM_OWNERSHIP_KEEP;
    if(pextract_face_vtx[i_part] == NULL) {
      owner = PDM_OWNERSHIP_KEEP;

      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      PDM_extract_part_connectivity_get(extrp,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                        &face_edge,
                                        &face_edge_idx,
                                        PDM_OWNERSHIP_KEEP);
      PDM_extract_part_connectivity_get(extrp,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx,
                                        &edge_vtx_idx,
                                        PDM_OWNERSHIP_KEEP);
      PDM_compute_face_vtx_from_face_and_edge(pn_extract_face[i_part], face_edge_idx, face_edge, edge_vtx, &pextract_face_vtx[i_part]);

      pextract_face_vtx_idx[i_part] = malloc((pn_extract_face[i_part]+1) * sizeof(int));
      for(int i = 0; i < pn_extract_face[i_part]+1; ++i) {
        pextract_face_vtx_idx[i_part][i] = face_edge_idx[i];
      }
    }


    PDM_extract_part_vtx_coord_get(extrp,
                                   i_part,
                                   &pextract_vtx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_CELL,
                                  &pextract_cell_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &pextract_face_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_VTX,
                                  &pextract_vtx_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    // log_trace(" %s --> %i \n", name_chr, pn_extract_cell[i_part]);
    /* Vtk en légende */
    char filename[999];
    sprintf(filename, "%s_%3.3d_%3.3d.vtk", name_chr, i_part, i_rank);
    PDM_vtk_write_polydata(filename,
                           pn_extract_vtx[i_part],
                           pextract_vtx[i_part],
                           pextract_vtx_ln_to_gn[i_part],
                           pn_extract_face[i_part],
                           pextract_face_vtx_idx[i_part],
                           pextract_face_vtx[i_part],
                           pextract_face_ln_to_gn[i_part],
                           NULL);

    if (owner == PDM_OWNERSHIP_KEEP) {
      free(pextract_face_vtx    [i_part]);
      free(pextract_face_vtx_idx[i_part]);
    }

  }

  // La visu concatene merge les valeurs donc on voit pas grand choses
  // _visu (name_chr,
  //        extrp->n_part_out,
  //        pn_extract_cell,
  //        pn_extract_face,
  //        pn_extract_vtx,
  //        pextract_cell_face_idx,
  //        pextract_cell_face,
  //        pextract_cell_ln_to_gn,
  //        pextract_face_vtx_idx,
  //        pextract_face_vtx,
  //        pextract_face_ln_to_gn,
  //        pextract_vtx,
  //        pextract_vtx_ln_to_gn);



  free(pn_extract_cell       );
  free(pn_extract_face       );
  free(pn_extract_vtx        );
  free(pextract_cell_face    );
  free(pextract_cell_face_idx);
  free(pextract_face_vtx     );
  free(pextract_face_vtx_idx );
  free(pextract_vtx          );
  free(pextract_cell_ln_to_gn);
  free(pextract_face_ln_to_gn);
  free(pextract_vtx_ln_to_gn );
}

static
void
_compute_extents_3d
(
       int      n_cell,
       int     *cell_face_idx,
       int     *cell_face,
       int     *face_vtx_idx,
       int     *face_vtx,
       double  *vtx_coord,
 const double   tolerance,
       double  *box_extents,
       double  *global_extents
)
{
  // const double tolerance   = 1.e-12;
  const double eps_extents = 1.e-7;
  const int dim = 3;

  /* Loop over cell */
  for(int i_cell = 0; i_cell < n_cell; ++i_cell ) {

    double *_extents = box_extents + 6 * i_cell;

    /* Init */
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[  i_dim] =  HUGE_VAL;
      _extents[3+i_dim] = -HUGE_VAL;
    }

    /* Loop over face and vtx to compute bbox */
    for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {

      int i_face = PDM_ABS(cell_face[idx_face])-1;

      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = face_vtx[idx_vtx]-1;

        for (int i_dim = 0; i_dim < 3; i_dim++) {
          double x = vtx_coord[3*i_vtx + i_dim];

          if (x < _extents[i_dim]) {
            _extents[i_dim] = x;
          }
          if (x > _extents[3+i_dim]) {
            _extents[3+i_dim] = x;
          }
        }
      }
    }

    double delta = 0.;
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      double x = _extents[3+i_dim] - _extents[i_dim];

      if (delta < x) {
        delta = x;
      }
    }

    if (delta > eps_extents) {
      delta *= tolerance;
    } else {
      delta = eps_extents;
    }

    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[i_dim]   -= delta;
      _extents[3+i_dim] += delta;
    }

    for (int k = 0; k < dim; k++) {
      global_extents[k]       = PDM_MIN(_extents[k    ], global_extents[k    ]);
      global_extents[dim + k] = PDM_MAX(_extents[dim+k], global_extents[dim+k]);
    }

  } /* End loop cell */
}


static
void
_compute_extents_2d_from_face_vtx
(
       int      n_face,
       int     *face_vtx_idx,
       int     *face_vtx,
       double  *vtx_coord,
 const double   tolerance,
       double  *box_extents,
       double  *global_extents
)
{
  // const double tolerance   = 1.e-12;
  const double eps_extents = 1.e-7;
  const int dim = 3;

  /* Loop over face */
  for(int i_face = 0; i_face < n_face; ++i_face ) {

    double *_extents = box_extents + 6 * i_face;

    /* Init */
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[  i_dim] =  HUGE_VAL;
      _extents[3+i_dim] = -HUGE_VAL;
    }

    /* Loop over face and vtx to compute bbox */
    for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
      int i_vtx = face_vtx[idx_vtx]-1;

      for (int i_dim = 0; i_dim < 3; i_dim++) {
        double x = vtx_coord[3*i_vtx + i_dim];

        if (x < _extents[i_dim]) {
          _extents[i_dim] = x;
        }
        if (x > _extents[3+i_dim]) {
          _extents[3+i_dim] = x;
        }
      }
    }

    double delta = 0.;
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      double x = _extents[3+i_dim] - _extents[i_dim];

      if (delta < x) {
        delta = x;
      }
    }

    if (delta > eps_extents) {
      delta *= tolerance;
    } else {
      delta = eps_extents;
    }

    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[i_dim]   -= delta;
      _extents[3+i_dim] += delta;
    }

    for (int k = 0; k < dim; k++) {
      global_extents[k]       = PDM_MIN(_extents[k    ], global_extents[k    ]);
      global_extents[dim + k] = PDM_MAX(_extents[dim+k], global_extents[dim+k]);
    }

  } /* End loop cell */
}


static
void
_compute_mesh_nodal_extents
(
        PDM_part_mesh_nodal_t   *mesh_nodal,
        int                      dim_mesh,
  const double                   tolerance,
        double                  *global_extents,
        double                ***extents_out
)
{
  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_MAX;
  switch (dim_mesh) {
  case 1:
    geom_kind = PDM_GEOMETRY_KIND_RIDGE;
    break;
  case 2:
    geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
    break;
  case 3:
    geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "invalid dimension %d\n", dim_mesh);
  }

  int n_part = PDM_part_mesh_nodal_n_part_get(mesh_nodal);

  int n_section = PDM_part_mesh_nodal_n_section_in_geom_kind_get(mesh_nodal, geom_kind);

  int *sections_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(mesh_nodal, geom_kind);

  double **extents = malloc(n_part * sizeof(double *));

  for (int i_part = 0; i_part < n_part; i_part++) {

    int part_n_elt = 0;

    int max_n_elt = 0;
    for (int isection = 0; isection < n_section; isection++) {
      int id_section_in_geom_kind = sections_id[isection];
      int id_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(mesh_nodal,
                                                                         geom_kind,
                                                                         id_section_in_geom_kind);
      int n_elt = PDM_part_mesh_nodal_section_n_elt_get(mesh_nodal,
                                                        id_section,
                                                        i_part);

      part_n_elt += n_elt;

      max_n_elt = PDM_MAX(max_n_elt, n_elt);
    }

    double *_extents = malloc(sizeof(double) * max_n_elt * 6);

    extents[i_part] = malloc(sizeof(double) * part_n_elt * 6);

    int idx = 0;
    for (int isection = 0; isection < n_section; isection++) {
      int id_section_in_geom_kind = sections_id[isection];
      int id_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(mesh_nodal,
                                                                         geom_kind,
                                                                         id_section_in_geom_kind);
      int *parent_num = PDM_part_mesh_nodal_section_parent_num_get(mesh_nodal,
                                                                   id_section,
                                                                   i_part,
                                                                   PDM_OWNERSHIP_KEEP);

      // PDM_g_num_t *_elt_g_num = PDM_part_mesh_nodal_section_g_num_get(mesh_nodal,
      //                                                               geom_kind,
      //                                                               id_section,
      //                                                               i_part);

      int n_elt = PDM_part_mesh_nodal_section_n_elt_get(mesh_nodal,
                                                        id_section,
                                                        i_part);

      PDM_part_mesh_nodal_section_elt_extents_compute(mesh_nodal,
                                                      id_section,
                                                      i_part,
                                                      tolerance,
                                                      _extents);

      for (int i = 0; i < n_elt; i++) {
        idx = i;
        if (parent_num) {
          idx = parent_num[i];
        }
        // part_elt_g_num[i_part][idx] = _elt_g_num[i];
        memcpy(extents[i_part] + 6*idx, _extents + 6*i, sizeof(double)*6);
      }
    }
    free(_extents);

    for (int i = 0; i < part_n_elt; i++) {
      for (int k = 0; k < 3; k++) {
        global_extents[k  ] = PDM_MIN(extents[i_part][6*i+k  ], global_extents[k  ]);
        global_extents[k+3] = PDM_MAX(extents[i_part][6*i+k+3], global_extents[k+3]);
      }
    }

  }


  *extents_out = extents;
}

static
void
_compute_part_mesh_extents
(
  PDM_part_mesh_t   *mesh,
  int                dim_mesh,
  const double       tolerance,
  double            *global_extents,
  double          ***extents_out
)
{
  int n_part = mesh->n_part;
  double **extents = malloc(n_part * sizeof(double *));
  if(dim_mesh == 3) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      int    *cell_face     = NULL;
      int    *cell_face_idx = NULL;
      int    *face_vtx_idx  = NULL;
      int    *face_vtx      = NULL;
      double *vtx_coord     = NULL;

      int n_cell = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL);

      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, &cell_face, &cell_face_idx, PDM_OWNERSHIP_BAD_VALUE);
      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx, PDM_OWNERSHIP_BAD_VALUE);

      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_BAD_VALUE); // Il faudrait un unchanged

      int *_face_vtx     = face_vtx;
      int *_face_vtx_idx = face_vtx_idx;
      extents[i_part] = malloc(6 * n_cell * sizeof(double));
      if (face_vtx == NULL) {
        int    *face_edge_idx  = NULL;
        int    *face_edge      = NULL;
        int    *edge_vtx_idx   = NULL;
        int    *edge_vtx       = NULL;
        int n_face = PDM_part_mesh_n_entity_get(mesh,
                                                i_part,
                                                PDM_MESH_ENTITY_FACE);
        PDM_part_mesh_connectivity_get(mesh,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                       &face_edge,
                                       &face_edge_idx,
                                       PDM_OWNERSHIP_BAD_VALUE);
        PDM_part_mesh_connectivity_get(mesh,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                       &edge_vtx,
                                       &edge_vtx_idx,
                                       PDM_OWNERSHIP_BAD_VALUE);
        PDM_compute_face_vtx_from_face_and_edge(n_face,
                                                face_edge_idx,
                                                face_edge,
                                                edge_vtx,
                                                &_face_vtx);
        _face_vtx_idx = face_edge_idx;
      }
      _compute_extents_3d(n_cell, cell_face_idx, cell_face,
                          _face_vtx_idx,
                          _face_vtx,
                          vtx_coord,
                          tolerance,
                          extents[i_part],
                          global_extents);
      if (face_vtx == NULL) {
        free(_face_vtx);
      }
    }
  } else if(dim_mesh == 2) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      int    *face_vtx      = NULL;
      int    *face_vtx_idx  = NULL;
      int    *face_edge_idx = NULL;
      int    *face_edge     = NULL;
      int    *edge_vtx_idx  = NULL;
      int    *edge_vtx      = NULL;
      double *vtx_coord     = NULL;

      // A gerer le cas mixte face_vtx ou face_edge + edge_vtx

      int n_face = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx, PDM_OWNERSHIP_BAD_VALUE);


      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_BAD_VALUE); // Il faudrait un unchanged
      extents[i_part] = malloc(6 * n_face * sizeof(double));

      if(face_vtx != NULL) {
        _compute_extents_2d_from_face_vtx(n_face,
                                          face_vtx_idx,
                                          face_vtx,
                                          vtx_coord,
                                          tolerance,
                                          extents[i_part],
                                          global_extents);
      } else {
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE , &face_edge, &face_edge_idx, PDM_OWNERSHIP_BAD_VALUE);
        assert(face_edge != NULL);
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx, &edge_vtx_idx, PDM_OWNERSHIP_BAD_VALUE);
        assert(edge_vtx_idx == NULL);
        int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);

        edge_vtx_idx = malloc((n_edge+1) * sizeof(int));
        for(int i_edge = 0; i_edge < n_edge+1; ++i_edge){
          edge_vtx_idx[i_edge] = 2 * i_edge;
        }
        _compute_extents_3d(n_face,
                            face_edge_idx,
                            face_edge,
                            edge_vtx_idx,
                            edge_vtx,
                            vtx_coord,
                            tolerance,
                            extents[i_part],
                            global_extents);
        free(edge_vtx_idx);
      }
    }

  } else {
    int    *edge_vtx_idx  = NULL;
    int    *edge_vtx      = NULL;
    double *vtx_coord     = NULL;
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);

      extents[i_part] = malloc(6 * n_edge * sizeof(double));
      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_BAD_VALUE); // Il faudrait un unchanged

      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx, &edge_vtx_idx, PDM_OWNERSHIP_BAD_VALUE);
      assert(edge_vtx_idx == NULL);
      edge_vtx_idx = malloc((n_edge+1) * sizeof(int));
      for(int i_edge = 0; i_edge < n_edge+1; ++i_edge){
        edge_vtx_idx[i_edge] = 2 * i_edge;
      }

      _compute_extents_2d_from_face_vtx(n_edge,
                                        edge_vtx_idx,
                                        edge_vtx,
                                        vtx_coord,
                                        tolerance,
                                        extents[i_part],
                                        global_extents);

      free(edge_vtx_idx);
    }
  }
  *extents_out = extents;
}

static
void
_select_elements_by_global_bbox
(
  PDM_part_mesh_t *mesh,
  int              dim_mesh,
  double         **box_extents,
  double          *g_mesh_global_extents,
  int            **n_extract_elmt_out,
  double        ***extract_box_extents_out,
  int           ***extract_elmt_init_location_out,
  PDM_g_num_t   ***extract_elmt_ln_to_gn_out
)
{
  int n_part = mesh->n_part;
  int i_rank;
  PDM_MPI_Comm_rank(mesh->comm, &i_rank);

  int          *n_extract_elmt             = malloc(n_part * sizeof(int         *));
  double      **extract_box_extents        = malloc(n_part * sizeof(double      *));
  int         **extract_elmt_init_location = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **extract_elmt_ln_to_gn      = malloc(n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_entity = 0;
    PDM_g_num_t* entity_ln_to_gn = NULL;
    if(dim_mesh == 3) {
      n_entity = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL);
      PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_CELL, &entity_ln_to_gn, PDM_OWNERSHIP_BAD_VALUE);
    } else if(dim_mesh == 2) {
      n_entity = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
      PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_FACE, &entity_ln_to_gn, PDM_OWNERSHIP_BAD_VALUE);
    } else if(dim_mesh == 1) {
      n_entity = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);
      PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_EDGE, &entity_ln_to_gn, PDM_OWNERSHIP_BAD_VALUE);
    }

    // char filename[999];
    // sprintf(filename, "titi_%i.vtk", n_entity);
    // PDM_vtk_write_boxes(filename,
    //                     n_entity,
    //                     box_extents[i_part],
    //                     NULL);

    n_extract_elmt[i_part] = 0;
    extract_box_extents       [i_part] = malloc(6 * n_entity * sizeof(double     ));
    extract_elmt_init_location[i_part] = malloc(3 * n_entity * sizeof(int        ));
    extract_elmt_ln_to_gn     [i_part] = malloc(    n_entity * sizeof(PDM_g_num_t));

    double *_box_extents = box_extents[i_part];

    for(int i = 0; i < n_entity; ++i) {

      double *box_min = _box_extents + 6*i;
      double *box_max = box_min + 3;

      int intersect = 1;
      for (int j = 0; j < 3; j++) {
        if (box_min[j] > g_mesh_global_extents[j+3] ||
            box_max[j] < g_mesh_global_extents[j  ]) {
          intersect = 0;
          break;
        }
      }

      if (intersect) {
        for (int j = 0; j < 6; j++) {
          extract_box_extents  [i_part][6*n_extract_elmt[i_part]+j] = _box_extents[6*i+j];
        }
        extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]  ] = i_rank;
        extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]+1] = i_part;
        extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]+2] = i;

        extract_elmt_ln_to_gn[i_part][n_extract_elmt[i_part]] = entity_ln_to_gn[i];

        n_extract_elmt[i_part]++;
      }
    }
    extract_box_extents       [i_part] = realloc(extract_box_extents       [i_part], 6 * n_entity * sizeof(double     ));
    extract_elmt_init_location[i_part] = realloc(extract_elmt_init_location[i_part], 3 * n_entity * sizeof(int        ));
    extract_elmt_ln_to_gn     [i_part] = realloc(extract_elmt_ln_to_gn     [i_part],     n_entity * sizeof(PDM_g_num_t));

  }

  *n_extract_elmt_out             = n_extract_elmt;
  *extract_box_extents_out        = extract_box_extents;
  *extract_elmt_init_location_out = extract_elmt_init_location;
  *extract_elmt_ln_to_gn_out      = extract_elmt_ln_to_gn;
}


static
void
_select_elements_by_global_bbox_nodal
(
       PDM_mesh_intersection_t   *mi,
 const int                        i_mesh,
       double                   **box_extents,
       double                    *g_mesh_global_extents,
       int                      **n_extract_elmt_out,
       double                  ***extract_box_extents_out,
       int                     ***extract_elmt_init_location_out,
       PDM_g_num_t             ***extract_elmt_ln_to_gn_out
 )
{
  PDM_MPI_Comm           comm       = mi->comm;
  PDM_part_mesh_nodal_t *mesh_nodal = mi->mesh_nodal[i_mesh];
  const int              dim_mesh   = mi->dim_mesh[i_mesh];

  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_MAX;
  switch (dim_mesh) {
  case 1:
    geom_kind = PDM_GEOMETRY_KIND_RIDGE;
    break;
  case 2:
    geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
    break;
  case 3:
    geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "invalid dimension %d\n", dim_mesh);
  }

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_part = PDM_part_mesh_nodal_n_part_get(mesh_nodal);

  int n_section = PDM_part_mesh_nodal_n_section_in_geom_kind_get(mesh_nodal, geom_kind);

  int *sections_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(mesh_nodal, geom_kind);

  int          *n_extract_elmt             = malloc(n_part * sizeof(int         *));
  double      **extract_box_extents        = malloc(n_part * sizeof(double      *));
  int         **extract_elmt_init_location = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **extract_elmt_ln_to_gn      = malloc(n_part * sizeof(PDM_g_num_t *));


  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_entity = 0;

    for (int isection = 0; isection < n_section; isection++) {
      int id_section_in_geom_kind = sections_id[isection];
      int id_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(mesh_nodal,
                                                                         geom_kind,
                                                                         id_section_in_geom_kind);

      int n_elt = PDM_part_mesh_nodal_section_n_elt_get(mesh_nodal,
                                                        id_section,
                                                        i_part);

      n_entity += n_elt;
    }

    n_extract_elmt[i_part] = 0;
    extract_box_extents       [i_part] = malloc(6 * n_entity * sizeof(double     ));
    extract_elmt_init_location[i_part] = malloc(3 * n_entity * sizeof(int        ));
    extract_elmt_ln_to_gn     [i_part] = malloc(    n_entity * sizeof(PDM_g_num_t));

    double *_box_extents = box_extents[i_part];

    for (int isection = 0; isection < n_section; isection++) {
      int id_section_in_geom_kind = sections_id[isection];
      int id_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(mesh_nodal,
                                                                         geom_kind,
                                                                         id_section_in_geom_kind);

      int n_elt = PDM_part_mesh_nodal_section_n_elt_get(mesh_nodal,
                                                        id_section,
                                                        i_part);

      int *parent_num = PDM_part_mesh_nodal_section_parent_num_get(mesh_nodal,
                                                                   id_section,
                                                                   i_part,
                                                                   PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *entity_ln_to_gn = PDM_part_mesh_nodal_g_num_get(mesh_nodal,
                                                                   id_section,
                                                                   i_part,
                                                                   PDM_OWNERSHIP_KEEP);

      for (int ielt = 0; ielt < n_elt; ielt++) {

        int i = ielt;
        if (parent_num != NULL) {
          i = parent_num[ielt];
        }

        double *box_min = _box_extents + 6*i;
        double *box_max = box_min + 3;

        int intersect = 1;
        for (int j = 0; j < 3; j++) {
          if (box_min[j] > g_mesh_global_extents[j+3] ||
              box_max[j] < g_mesh_global_extents[j  ]) {
            intersect = 0;
            break;
          }
        }

        if (intersect) {
          for (int j = 0; j < 6; j++) {
            extract_box_extents  [i_part][6*n_extract_elmt[i_part]+j] = _box_extents[6*i+j];
          }
          extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]  ] = i_rank;
          extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]+1] = i_part;
          extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]+2] = i;

          extract_elmt_ln_to_gn[i_part][n_extract_elmt[i_part]] = entity_ln_to_gn[ielt];

          n_extract_elmt[i_part]++;
        }
      }

    }

    extract_box_extents       [i_part] = realloc(extract_box_extents       [i_part], 6 * n_entity * sizeof(double     ));
    extract_elmt_init_location[i_part] = realloc(extract_elmt_init_location[i_part], 3 * n_entity * sizeof(int        ));
    extract_elmt_ln_to_gn     [i_part] = realloc(extract_elmt_ln_to_gn     [i_part],     n_entity * sizeof(PDM_g_num_t));
  }

  *n_extract_elmt_out             = n_extract_elmt;
  *extract_box_extents_out        = extract_box_extents;
  *extract_elmt_init_location_out = extract_elmt_init_location;
  *extract_elmt_ln_to_gn_out      = extract_elmt_ln_to_gn;
}

static
void
_redistrib_boxes
(
 PDM_MPI_Comm    comm,
 PDM_box_set_t  *boxes_mesh_a,
 PDM_box_set_t  *boxes_mesh_b,
 int            *box_a_to_box_b_idx,
 int            *box_a_to_box_b,
 int           **redistribute_box_a_to_box_b_idx,
 int           **redistribute_box_a_to_box_b
)
{

  int n_elt_mesh_a = PDM_box_set_get_size (boxes_mesh_a);
  int n_elt_mesh_b = PDM_box_set_get_size (boxes_mesh_b);

  PDM_g_num_t *gnum_elt_mesh_a = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_mesh_a);
  PDM_g_num_t *gnum_elt_mesh_b = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_mesh_b);

  /*****************************************************************************
   *                                                                           *
   *  Transfer intersection information from partitions to blocks              *
   * with PDM_part_to_block_exch function                                      *
   *                                                                           *
   *  Results :                                                                *
   *      - block_a_boxes_b_idx                                                *
   *      - block_a_boxes_b_gnum_data                                          *
   *                                                                           *
   ****************************************************************************/

  /*
   * Tentative Bruno :
   *   - Ponderate work
   *   - Extract only cell with job
   */
  double* weight = (double *) malloc( n_elt_mesh_a * sizeof(double));
  // PDM_g_num_t* extract_mesh_a_g_num = (PDM_g_num_t *) malloc( n_elt_mesh_a * sizeof(PDM_g_num_t));
  for (int i = 0; i < n_elt_mesh_a; i++) {
    weight[i] = box_a_to_box_b_idx[i+1] - box_a_to_box_b_idx[i];
  }

  // TODO : Geometric to better locality
  PDM_part_to_block_t *ptb_boxes_a = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_MERGE,
                                                              1.,
                                            (PDM_g_num_t **) &gnum_elt_mesh_a,
                                                             &weight,
                                                              &n_elt_mesh_a,
                                                              1,
                                                              comm);

  int n_elt_block_a = PDM_part_to_block_n_elt_block_get (ptb_boxes_a);
  free(weight);

  PDM_g_num_t *block_gnum_a = PDM_part_to_block_block_gnum_get (ptb_boxes_a);

  int *part_stride_a = (int *) malloc (sizeof(int) * n_elt_mesh_a);

  for (int i = 0; i < n_elt_mesh_a; i++) {
    part_stride_a[i] = box_a_to_box_b_idx[i+1] - box_a_to_box_b_idx[i];
  }

  /*
   * Exchange connectivity box_a_to_box_b
   */
  PDM_g_num_t *box_a_to_box_b_g_num = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) *
                                                              box_a_to_box_b_idx[n_elt_mesh_a]);

  for (int k = 0; k < box_a_to_box_b_idx[n_elt_mesh_a]; k++) {
    box_a_to_box_b_g_num[k] = gnum_elt_mesh_b[box_a_to_box_b[k]];
  }

  int         *block_a_boxes_b_stride;
  PDM_g_num_t *block_a_boxes_b_gnum_data;

  PDM_part_to_block_exch (ptb_boxes_a,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &part_stride_a,
               (void **) &box_a_to_box_b_g_num,
                         &block_a_boxes_b_stride,
               (void **) &block_a_boxes_b_gnum_data);
  free(box_a_to_box_b_g_num);


  // if (1) {
  //   int idx = 0;
  //   log_trace("--- Avant ---\n");
  //   for (int i = 0; i < n_elt_block_a; i++) {
  //     log_trace("faceA "PDM_FMT_G_NUM": faceB ", block_gnum_a[i]);
  //     for (int j = 0; j < block_a_boxes_b_stride[i]; j++) {
  //       log_trace(" "PDM_FMT_G_NUM, block_a_boxes_b_gnum_data[idx++]);
  //     }
  //     log_trace("\n");
  //   }
  // }

  /* Remove duplicates */
  if (1) {
    int idx_read  = 0;
    int idx_write = 0;
    for (int i = 0; i < n_elt_block_a; i++) {

      int n = block_a_boxes_b_stride[i];

      int m = PDM_inplace_unique_long(block_a_boxes_b_gnum_data + idx_read,
                                      NULL,
                                      0,
                                      n-1);

      for (int j = 0; j < m; j++) {
        block_a_boxes_b_gnum_data[idx_write++] = block_a_boxes_b_gnum_data[idx_read+j];
      }
      block_a_boxes_b_stride[i] = m;

      idx_read += n;
    }
  }

  // if (1) {
  //   int idx = 0;
  //   log_trace("--- Après ---\n");
  //   for (int i = 0; i < n_elt_block_a; i++) {
  //     log_trace("faceA "PDM_FMT_G_NUM": faceB ", block_gnum_a[i]);
  //     for (int j = 0; j < block_a_boxes_b_stride[i]; j++) {
  //       log_trace(" "PDM_FMT_G_NUM, block_a_boxes_b_gnum_data[idx++]);
  //     }
  //     log_trace("\n");
  //   }
  // }




  /*****************************************************************************
   *                                                                           *
   * Redistribute boxes_mesh_a intersections to ensure a good load balacing          *
   * in comm MPI communicator                                              *
   * This step removes intersections found many times on different ranks       *
   *                                                                           *
   * After this step, data are stored in a block with n_elt_block_a, block_gnum_a,
   * part_stride_a                                                              *
   *                                                                           *
   ****************************************************************************/
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  /*
   * Redistribute data boxes_mesh_a from blockB distribution with a PDM_box_distrib_t
   * structure
   * TODO: - Build a new PDM_box_distrib_t structure more simple
   *       - Hide PDM_box_distrib_t attributes
   */

  PDM_l_num_t *destination = PDM_part_to_block_destination_get (ptb_boxes_a);

  PDM_g_num_t n_g_elmt_mesh_a = PDM_box_set_get_global_size(boxes_mesh_a);

  PDM_box_distrib_t *distrib_a = PDM_box_distrib_create(n_elt_mesh_a,
                                                        n_g_elmt_mesh_a,
                                                        1, // Don't use in this case
                                                        comm);

  PDM_g_num_t n_g_elmt_mesh_b = PDM_box_set_get_global_size(boxes_mesh_b);

  PDM_box_distrib_t *distrib_b = PDM_box_distrib_create(n_elt_mesh_b,
                                                        n_g_elmt_mesh_b,
                                                        1, // Don't use in this case
                                                        comm);


  int *count_elts_a = (int *) malloc (sizeof(int) * n_rank);
  int *count_elts_b = (int *) malloc (sizeof(int) * n_rank);

  for (int i = 0; i < n_rank + 1; i++) {
    distrib_a->index[i] = 0;
    distrib_b->index[i] = 0;
  }

  for (int i = 0; i < n_rank; i++) {
    count_elts_a[i] = 0;
    count_elts_b[i] = 0;
  }

  for (int i = 0; i < n_elt_mesh_a; i++) {
    int t_rank = destination[i] + 1;

    if(part_stride_a[i] > 0 ) { // To see with Eric and Bastien --> I use it to extract only the intersect part
      distrib_a->index[t_rank]++;
      distrib_b->index[t_rank] += part_stride_a[i];
    }
  }

  for (int i = 0; i < n_rank; i++) {
    distrib_a->index[i+1] += distrib_a->index[i];
    distrib_b->index[i+1] += distrib_b->index[i];
  }

  distrib_a->list = (int *) malloc (sizeof(int) * distrib_a->index[n_rank]);
  distrib_b->list = (int *) malloc (sizeof(int) * distrib_b->index[n_rank]);

  for (int i = 0; i < n_elt_mesh_a; i++) {
    if(part_stride_a[i] > 0 ) { // To see with Eric and Bastien --> I use it to extract only the intersect part
      int t_rank = destination[i]; // EQU + 1; mais ce n est pas necessaire
      int idx_a = distrib_a->index[t_rank] + (count_elts_a[t_rank]++);
      distrib_a->list[idx_a] = i;
      int idx_b = distrib_b->index[t_rank] + count_elts_b[t_rank];
      count_elts_b[t_rank] += part_stride_a[i];
      int k=0;
      for (int j = box_a_to_box_b_idx[i]; j < box_a_to_box_b_idx[i+1]; j++) {
        distrib_b->list[idx_b+k++] = box_a_to_box_b[j];
      }
    }
  }


  free (part_stride_a);
  free (count_elts_a);
  free (count_elts_b);

  PDM_box_distrib_clean (distrib_a);
  PDM_box_distrib_clean (distrib_b);

  PDM_box_set_redistribute (distrib_a, boxes_mesh_a);
  PDM_box_set_redistribute (distrib_b, boxes_mesh_b);

  PDM_box_distrib_destroy (&distrib_a);
  PDM_box_distrib_destroy (&distrib_b);

  PDM_box_set_remove_duplicate (boxes_mesh_a);
  PDM_box_set_remove_duplicate (boxes_mesh_b);

  /*
   * All boxes are redistribute we need to update box_a_to_box_b array
   *    - Caution if morton / hilbert the array block_gnum_a IS NOT order
   */
  n_elt_mesh_a    = PDM_box_set_get_size (boxes_mesh_a); // Caution not the same of the fist call because redistibute
  n_elt_mesh_b    = PDM_box_set_get_size (boxes_mesh_b); // Caution not the same of the fist call because redistibute

  gnum_elt_mesh_a = (PDM_g_num_t *) PDM_box_set_get_g_num(boxes_mesh_a);
  gnum_elt_mesh_b = (PDM_g_num_t *) PDM_box_set_get_g_num(boxes_mesh_b);

  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(block_gnum_a,
                                                                        n_elt_block_a,
                                                (const PDM_g_num_t **)  &gnum_elt_mesh_a,
                                                                        &n_elt_mesh_a,
                                                                        1,
                                                                        comm);

  PDM_part_to_block_free (ptb_boxes_a);
  int         **tmp_redistribute_box_a_to_box_b_n     = NULL;
  PDM_g_num_t **tmp_redistribute_box_a_to_box_b_g_num = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         block_a_boxes_b_stride,
                         block_a_boxes_b_gnum_data,
                         &tmp_redistribute_box_a_to_box_b_n,
           (void ***)    &tmp_redistribute_box_a_to_box_b_g_num);
  free (block_a_boxes_b_stride);
  free (block_a_boxes_b_gnum_data);

  // if (1) {
  //   int idx = 0;
  //   log_trace("--- 2 ---\n");
  //   for (int i = 0; i < n_elt_mesh_a; i++) {
  //     log_trace("faceA "PDM_FMT_G_NUM": faceB ", gnum_elt_mesh_a[i]);
  //     for (int j = 0; j < tmp_redistribute_box_a_to_box_b_n[0][i]; j++) {
  //       log_trace(" "PDM_FMT_G_NUM, tmp_redistribute_box_a_to_box_b_g_num[0][idx++]);
  //     }
  //     log_trace("\n");
  //   }
  // }

  PDM_block_to_part_free(btp);

  int         *redistribute_box_a_to_box_b_n     = tmp_redistribute_box_a_to_box_b_n    [0];
  PDM_g_num_t *redistribute_box_a_to_box_b_g_num = tmp_redistribute_box_a_to_box_b_g_num[0];
  free(tmp_redistribute_box_a_to_box_b_n    );
  free(tmp_redistribute_box_a_to_box_b_g_num);


  /*
   * Translate in frame of B
   */

  int         *order              = (int         *) malloc(n_elt_mesh_b * sizeof(int        ));
  PDM_g_num_t *gnum_elt_mesh_b_cp = (PDM_g_num_t *) malloc(n_elt_mesh_b * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_elt_mesh_b; ++i ) {
    order             [i] = i;
    gnum_elt_mesh_b_cp[i] = gnum_elt_mesh_b[i];
  }

  PDM_sort_long(gnum_elt_mesh_b_cp, order, n_elt_mesh_b);


  int *_redistribute_box_a_to_box_b_idx = (int *) malloc((n_elt_mesh_a+1) * sizeof(int));
  _redistribute_box_a_to_box_b_idx[0] = 0;
  int n_tot_connect = 0;
  for(int i = 0; i < n_elt_mesh_a; ++i) {
    n_tot_connect += redistribute_box_a_to_box_b_n[i];
    _redistribute_box_a_to_box_b_idx[i+1] = _redistribute_box_a_to_box_b_idx[i] + redistribute_box_a_to_box_b_n[i];
  }

  int *_redistribute_box_a_to_box_b = (int *) malloc( n_tot_connect * sizeof(int));

  for(int i = 0; i < n_tot_connect; ++i) {
    int pos = PDM_binary_search_long(redistribute_box_a_to_box_b_g_num[i], gnum_elt_mesh_b_cp, n_elt_mesh_b);
    _redistribute_box_a_to_box_b[i] = order[pos];
  }

  *redistribute_box_a_to_box_b_idx = _redistribute_box_a_to_box_b_idx;
  *redistribute_box_a_to_box_b     = _redistribute_box_a_to_box_b;

  free(order);
  free(gnum_elt_mesh_b_cp);
  free(redistribute_box_a_to_box_b_n    );
  free(redistribute_box_a_to_box_b_g_num);
}


static
PDM_extract_part_t*
_create_extract_part
(
  PDM_part_mesh_t              *mesh,
  int                           dim_mesh,
  PDM_mesh_intersection_kind_t  intersect_kind,
  PDM_box_set_t                *boxes_meshes
)
{
  PDM_g_num_t *gnum_elt_mesh = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_meshes);
  int n_elt_mesh = PDM_box_set_get_size (boxes_meshes);

  int n_part_out = 1;
  PDM_extract_part_t* extrp_mesh = PDM_extract_part_create(dim_mesh,
                                                           mesh->n_part,
                                                           n_part_out,
                                                           PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                                           PDM_SPLIT_DUAL_WITH_HILBERT, // Not used
                                                           PDM_FALSE,                   // compute_child_gnum
                                                           PDM_OWNERSHIP_KEEP,
                                                           mesh->comm);

  PDM_g_num_t *target_g_num = gnum_elt_mesh;
  if (intersect_kind == PDM_MESH_INTERSECTION_KIND_PREPROCESS) { 
    target_g_num = malloc(sizeof(PDM_g_num_t) * n_elt_mesh);
    memcpy(target_g_num, gnum_elt_mesh, sizeof(PDM_g_num_t) * n_elt_mesh);
    PDM_extract_part_target_gnum_keep_ownnership(extrp_mesh);
  }


  // printf("n_elt_mesh = %i  \n", n_elt_mesh);


  int *init_location_elt_mesh = (int  *) PDM_box_set_origin_get(boxes_meshes);


  for(int i_part = 0; i_part < mesh->n_part; ++i_part) {

    int n_cell = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL);
    int n_face = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
    int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);
    int n_vtx  = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_VTX);

    PDM_g_num_t *cell_ln_to_gn = NULL;
    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *edge_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;

    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_CELL, &cell_ln_to_gn, PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_FACE, &face_ln_to_gn, PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_EDGE, &edge_ln_to_gn, PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_VTX,  &vtx_ln_to_gn , PDM_OWNERSHIP_BAD_VALUE);

    int *cell_face     = NULL;
    int *cell_face_idx = NULL;
    int *face_vtx      = NULL;
    int *face_vtx_idx  = NULL;
    int *face_edge     = NULL;
    int *face_edge_idx = NULL;
    int *edge_vtx      = NULL;
    int *edge_vtx_idx  = NULL;

    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, &cell_face, &cell_face_idx, PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx , PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE, &face_edge, &face_edge_idx, PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , &edge_vtx , &edge_vtx_idx , PDM_OWNERSHIP_BAD_VALUE);

    double *vtx_coord = NULL;
    PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_BAD_VALUE);

    PDM_extract_part_part_set(extrp_mesh,
                              i_part,
                              n_cell,
                              n_face,
                              n_edge,
                              n_vtx ,
                              cell_face_idx,
                              cell_face,
                              face_edge_idx,
                              face_edge,
                              edge_vtx,
                              face_vtx_idx,
                              face_vtx,
                              cell_ln_to_gn,
                              face_ln_to_gn,
                              edge_ln_to_gn,
                              vtx_ln_to_gn,
                              vtx_coord);
  }




  /*  Setup target frame */
  PDM_extract_part_target_set(extrp_mesh, 0, n_elt_mesh, target_g_num, init_location_elt_mesh);

  PDM_extract_part_compute(extrp_mesh);

  return extrp_mesh;
}


static
PDM_extract_part_t*
_create_extract_part_nodal
(
 PDM_mesh_intersection_t      *mi,
 const int                     i_mesh,
 PDM_box_set_t                 *boxes_meshes
)
{
  PDM_MPI_Comm comm = mi->comm;
  PDM_part_mesh_nodal_t *mesh_nodal = mi->mesh_nodal[i_mesh];
  const int dim_mesh = mi->dim_mesh[i_mesh];

  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_MAX;
  switch (dim_mesh) {
  case 1:
    geom_kind = PDM_GEOMETRY_KIND_RIDGE;
    break;
  case 2:
    geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
    break;
  case 3:
    geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "invalid dimension %d\n", dim_mesh);
  }

  int n_part       = mi->n_part_mesh[i_mesh];
  int n_section    = 0;
  int *sections_id = NULL;
  PDM_part_mesh_nodal_elmts_t *pmne = NULL;
  if (mesh_nodal != NULL) {
    switch (dim_mesh) {
    case 1:
      pmne = mesh_nodal->ridge;
      break;
    case 2:
      pmne = mesh_nodal->surfacic;
      break;
    case 3:
      pmne = mesh_nodal->volumic;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid dimension %d\n", dim_mesh);
    }

    assert(n_part == PDM_part_mesh_nodal_n_part_get(mesh_nodal));
    n_section   = PDM_part_mesh_nodal_n_section_in_geom_kind_get  (mesh_nodal, geom_kind);
    sections_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(mesh_nodal, geom_kind);
  }

  PDM_part_mesh_nodal_elmts_extend_to_encompassing_comm(comm,
                                                        n_part,
                                                        &pmne);

  PDM_g_num_t *gnum_elt_mesh = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_meshes);

  int n_elt_mesh = PDM_box_set_get_size (boxes_meshes);


  int n_part_out = 1;
  PDM_extract_part_t* extrp_mesh = PDM_extract_part_create(dim_mesh,
                                                           n_part,
                                                           n_part_out,
                                                           PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                                           PDM_SPLIT_DUAL_WITH_HILBERT, // Not used
                                                           PDM_TRUE,//PDM_FALSE,                   // compute_child_gnum
                                                           PDM_OWNERSHIP_KEEP,
                                                           comm);

  PDM_g_num_t *target_g_num = gnum_elt_mesh;
  if (mi->intersect_kind == PDM_MESH_INTERSECTION_KIND_PREPROCESS) { 
    target_g_num = malloc(sizeof(PDM_g_num_t) * n_elt_mesh);
    memcpy(target_g_num, gnum_elt_mesh, sizeof(PDM_g_num_t) * n_elt_mesh);
    PDM_extract_part_target_gnum_keep_ownnership(extrp_mesh);
  }

  // printf("n_elt_mesh = %i  \n", n_elt_mesh);


  int *init_location_elt_mesh = (int  *) PDM_box_set_origin_get(boxes_meshes);


  PDM_extract_part_part_nodal_set(extrp_mesh, pmne);

  /* Set vtx coordinates */
  for (int i_part = 0; i_part < n_part; ++i_part) {
    int n_vtx = 0;
    double      *vtx_coord    = NULL;
    PDM_g_num_t *vtx_ln_to_gn = NULL;

    if (mesh_nodal != NULL) {
      n_vtx = PDM_part_mesh_nodal_n_vtx_get(mesh_nodal, i_part);
      vtx_coord    = (double      *) PDM_part_mesh_nodal_vtx_coord_get(mesh_nodal, i_part);
      vtx_ln_to_gn = (PDM_g_num_t *) PDM_part_mesh_nodal_vtx_g_num_get(mesh_nodal, i_part);
    }

    int part_n_elt = 0;
    for (int isection = 0; isection < n_section; isection++) {
      int id_section_in_geom_kind = sections_id[isection];
      int id_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(mesh_nodal,
                                                                         geom_kind,
                                                                         id_section_in_geom_kind);

      part_n_elt += PDM_part_mesh_nodal_section_n_elt_get(mesh_nodal,
                                                          id_section,
                                                          i_part);
    }


    int n_cell = 0;
    int n_face = 0;
    int n_edge = 0;
    PDM_g_num_t *cell_ln_to_gn = NULL;
    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *edge_ln_to_gn = NULL;

    switch (dim_mesh) {
    case 1:
      n_edge = part_n_elt;
      break;
    case 2:
      n_face = part_n_elt;
      break;
    case 3:
      n_cell = part_n_elt;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "invalid dimension %d\n", dim_mesh);
    }


    PDM_extract_part_part_set(extrp_mesh,
                              i_part,
                              n_cell,
                              n_face,
                              n_edge,
                              n_vtx,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              cell_ln_to_gn,
                              face_ln_to_gn,
                              edge_ln_to_gn,
                              vtx_ln_to_gn,
                              vtx_coord);
  }

  /*  Setup target frame */
  PDM_extract_part_target_set(extrp_mesh, 0, n_elt_mesh, target_g_num, init_location_elt_mesh);
  // PDM_g_num_t *target_g_num = malloc(sizeof(PDM_g_num_t) * n_elt_mesh);
  // memcpy(target_g_num, gnum_elt_mesh, sizeof(PDM_g_num_t) * n_elt_mesh);
  // PDM_extract_part_target_set(extrp_mesh, 0, n_elt_mesh, target_g_num, init_location_elt_mesh);

  PDM_extract_part_compute(extrp_mesh);

  return extrp_mesh;
}



 static PDM_ownership_t
_get_extracted_mesh_vol
(
 PDM_extract_part_t  *extrp,
 int                 *n_cell,
 int                 *n_face,
 int                 *n_vtx,
 int                **cell_face_idx,
 int                **cell_face,
 int                **face_vtx_idx,
 int                **face_vtx,
 double             **vtx_coord,
 PDM_g_num_t        **cell_ln_to_gn,
 PDM_g_num_t        **vtx_ln_to_gn
 )
 {
  *n_cell = PDM_extract_part_connectivity_get(extrp,
                                              0,
                                              PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                              cell_face,
                                              cell_face_idx,
                                              PDM_OWNERSHIP_KEEP);

  *n_face = PDM_extract_part_connectivity_get(extrp,
                                              0,
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                              face_vtx,
                                              face_vtx_idx,
                                              PDM_OWNERSHIP_KEEP);

  PDM_ownership_t owner = PDM_OWNERSHIP_KEEP;
  if (*face_vtx == NULL) {
    owner = PDM_OWNERSHIP_USER;

    int *face_edge = NULL;
    PDM_extract_part_connectivity_get(extrp,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                      &face_edge,
                                      face_vtx_idx,
                                      PDM_OWNERSHIP_KEEP);

    int *edge_vtx     = NULL;
    int *edge_vtx_idx = NULL;
    PDM_extract_part_connectivity_get(extrp,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &edge_vtx,
                                      &edge_vtx_idx,
                                      PDM_OWNERSHIP_KEEP);

    PDM_compute_face_vtx_from_face_and_edge(*n_face,
                                            *face_vtx_idx,
                                            face_edge,
                                            edge_vtx,
                                            face_vtx);
  }

  *n_vtx = PDM_extract_part_vtx_coord_get(extrp,
                                          0,
                                          vtx_coord,
                                          PDM_OWNERSHIP_KEEP);

  *cell_ln_to_gn = extrp->target_gnum[0];
  // PDM_extract_part_parent_ln_to_gn_get(extrp,
  //                                      0,
  //                                      PDM_MESH_ENTITY_CELL,
  //                                      cell_ln_to_gn,
  //                                      PDM_OWNERSHIP_KEEP);

  PDM_extract_part_parent_ln_to_gn_get(extrp,
                                       0,
                                       PDM_MESH_ENTITY_VTX,
                                       vtx_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);

  return owner;
 }


static void
_export_ensight3d
(
 PDM_MPI_Comm  comm,
 const char   *name,
 int           n_cell,
 int           n_face,
 int           n_vtx,
 int          *cell_face_idx,
 int          *cell_face,
 int          *face_vtx_idx,
 int          *face_vtx,
 double       *vtx_coord,
 PDM_g_num_t  *cell_ln_to_gn,
 PDM_g_num_t  *vtx_ln_to_gn
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_gen_gnum_t *gnum_vtx = PDM_gnum_create(3, 1, PDM_FALSE, 1., comm, PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_from_parents(gnum_vtx,
                            0,
                            n_vtx,
                            vtx_ln_to_gn);

  PDM_gnum_compute(gnum_vtx);

  PDM_g_num_t *extract_vtx_ln_to_gn = PDM_gnum_get(gnum_vtx,
                                                   0);


  PDM_gen_gnum_t *gnum_cell = PDM_gnum_create(3, 1, PDM_FALSE, 1., comm, PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_from_parents(gnum_cell,
                            0,
                            n_cell,
                            cell_ln_to_gn);

  PDM_gnum_compute(gnum_cell);

  PDM_g_num_t *extract_cell_ln_to_gn = PDM_gnum_get(gnum_cell,
                                                    0);


  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        name,
                                        name,
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       "iso_surface_volume_mesh",
                                       1);

  int id_var_rank = PDM_writer_var_create(wrt,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "i_rank");

  int id_var_gnum = PDM_writer_var_create(wrt,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "cell_g_num");

  int id_var_vol = PDM_writer_var_create(wrt,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "cell_volume");


  PDM_writer_step_beg(wrt, 0.);

  // int *cell_face_n = malloc(sizeof(int) * n_cell);
  // for (int i = 0; i < n_cell; i++) {
  //   cell_face_n[i] = cell_face_idx[i+1] - cell_face_idx[i];
  // }

  // int *face_vtx_n  = malloc(sizeof(int) * n_face);
  // for (int i = 0; i < n_face; i++) {
  //   face_vtx_n[i] = face_vtx_idx[i+1] - face_vtx_idx[i];
  // }


  PDM_writer_geom_coord_set(wrt,
                            id_geom,
                            0,
                            n_vtx,
                            vtx_coord,
                            extract_vtx_ln_to_gn,
                            PDM_OWNERSHIP_USER);

  PDM_writer_geom_cell3d_cellface_add(wrt,
                                      id_geom,
                                      0,
                                      n_cell,
                                      n_face,
                                      face_vtx_idx,
                                      NULL,//face_vtx_n,
                                      face_vtx,
                                      cell_face_idx,
                                      NULL,//cell_face_n,
                                      cell_face,
                                      extract_cell_ln_to_gn);

  PDM_writer_geom_write(wrt,
                        id_geom);

  PDM_real_t *val_rank = malloc(sizeof(PDM_real_t) * n_cell);
  PDM_real_t *val_gnum = malloc(sizeof(PDM_real_t) * n_cell);
  PDM_real_t *val_vol  = malloc(sizeof(PDM_real_t) * n_cell);

  double *volume = malloc(sizeof(double) * n_cell);
  double *center = malloc(sizeof(double) * n_cell * 3);
  PDM_geom_elem_polyhedra_properties_triangulated(1,
                                                  n_cell,
                                                  n_face,
                                                  face_vtx_idx,
                                                  face_vtx,
                                                  cell_face_idx,
                                                  cell_face,
                                                  n_vtx,
                                                  vtx_coord,
                                                  volume,
                                                  center,
                                                  NULL,
                                                  NULL);

  for (int i = 0; i < n_cell; i++) {
    val_rank[i] = (PDM_real_t) i_rank;
    val_gnum[i] = (PDM_real_t) cell_ln_to_gn[i];
    val_vol [i] = (PDM_real_t) volume[i];
  }
  free(volume);
  free(center);

  PDM_writer_var_set(wrt,
                     id_var_rank,
                     id_geom,
                     0,
                     val_rank);

  PDM_writer_var_write(wrt,
                       id_var_rank);
  PDM_writer_var_free(wrt,
                      id_var_rank);


  PDM_writer_var_set(wrt,
                     id_var_gnum,
                     id_geom,
                     0,
                     val_gnum);

  PDM_writer_var_write(wrt,
                       id_var_gnum);
  PDM_writer_var_free(wrt,
                      id_var_gnum);

  PDM_writer_var_set(wrt,
                     id_var_vol,
                     id_geom,
                     0,
                     val_vol);

  PDM_writer_var_write(wrt,
                       id_var_vol);
  PDM_writer_var_free(wrt,
                      id_var_vol);

  PDM_writer_step_end(wrt);

  PDM_writer_free(wrt);

  free(val_rank);
  free(val_gnum);
  free(val_vol );
  // free(cell_face_n);
  // free(face_vtx_n);

  PDM_gnum_free(gnum_vtx );
  PDM_gnum_free(gnum_cell);
}


static void
_build_ptp
(
 PDM_mesh_intersection_t *mi,
 int                     *elt_a_elt_b_idx,
 int                     *elt_a_elt_b,
 double                  *elt_a_elt_b_volume
 )
{
  PDM_mpi_comm_kind_t comm_kind = PDM_MPI_COMM_KIND_P2P;//COLLECTIVE;

  PDM_extract_part_t *extrp_mesh_a = mi->extrp_mesh[0];
  PDM_extract_part_t *extrp_mesh_b = mi->extrp_mesh[1];


  PDM_mesh_entities_t entity_type[2];
  for (int imesh = 0; imesh < 2; imesh++) {
    if (mi->dim_mesh[imesh] == 1) {
      entity_type[imesh] = PDM_MESH_ENTITY_EDGE;
    }
    else if (mi->dim_mesh[imesh] == 2) {
      entity_type[imesh] = PDM_MESH_ENTITY_FACE;
    }
    else if (mi->dim_mesh[imesh] == 3) {
      entity_type[imesh] = PDM_MESH_ENTITY_CELL;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "invalid dimension for mesh %d\n", imesh);
    }
  }


  // PDM_g_num_t *elt_a_ln_to_gn = NULL;
  int n_elt_a = extrp_mesh_a->n_target[0];
  // elt_a_ln_to_gn = extrp_mesh_a->target_gnum[0];
  // int n_elt_a = PDM_extract_part_parent_ln_to_gn_get(extrp_mesh_a,
  //                                                    0,
  //                                                    entity_type_a,
  //                                                    &elt_a_ln_to_gn,
  //                                                    PDM_OWNERSHIP_KEEP);


  PDM_g_num_t *elt_b_ln_to_gn = NULL;
  int n_elt_b = extrp_mesh_b->n_target[0];
  elt_b_ln_to_gn = extrp_mesh_b->target_gnum[0];
  // it n_elt_b = PDM_extractpart_parent_ln_to_gn_get(extrp_mesh_b,
  //                                                    0,
  //                                                    entity_type_b,
  //                                                    &elt_b_ln_to_gn,
  //                                                    PDM_OWNERSHIP_KEEP);

  /* Get all init locations of extracted faces B */
  // may not work if multiple init locations...

  int *elt_b_init_loc_n = PDM_array_const_int(n_elt_b, 1);
  int *elt_b_init_loc   = extrp_mesh_b->target_location[0];
  int *elt_b_init_loc_idx = PDM_array_new_idx_from_sizes_int(elt_b_init_loc_n, n_elt_b);

  /* Remove false positives */
  // !!! do not modify elt_a_elt_b(_idx) (owned by someone else)
  int idx_read  = 0;
  int idx_write = 0;
  int s_elt_a_elt_b_init_loc = elt_a_elt_b_idx[n_elt_a] * 2;
  int idx_write_init_loc = 0;
  int         *elt_a_elt_b_n          = malloc(sizeof(int        ) * n_elt_a);
  PDM_g_num_t *elt_a_elt_b_g_num      = malloc(sizeof(PDM_g_num_t) * elt_a_elt_b_idx[n_elt_a]);
  int         *elt_a_elt_b_init_loc_n = malloc(sizeof(int        ) * elt_a_elt_b_idx[n_elt_a]);
  int         *elt_a_elt_b_init_loc   = malloc(sizeof(int        ) * s_elt_a_elt_b_init_loc * 3);
  int         *elt_a_elt_b_init_loc_stride = malloc(sizeof(int) * n_elt_a);

  for (int elt_a_id = 0; elt_a_id < n_elt_a; elt_a_id++) {
    int n = elt_a_elt_b_idx[elt_a_id+1] - idx_read;
    elt_a_elt_b_n[elt_a_id] = 0;
    elt_a_elt_b_init_loc_stride[elt_a_id] = 0;
    for (int i = 0; i < n; i++) {
      if (elt_a_elt_b_volume[idx_read+i] > 1e-15) {
        // elt_a_elt_b       [idx_write] = elt_a_elt_b       [idx_read+i];
        int elt_b_id = elt_a_elt_b[idx_read+i];
        elt_a_elt_b_g_num [idx_write]     = elt_b_ln_to_gn[elt_b_id];
        elt_a_elt_b_volume[idx_write]     = elt_a_elt_b_volume[idx_read+i];
        elt_a_elt_b_init_loc_n[idx_write] = elt_b_init_loc_n[elt_b_id];
        if (idx_write_init_loc + elt_b_init_loc[elt_b_id] >= s_elt_a_elt_b_init_loc) {
          s_elt_a_elt_b_init_loc = PDM_MAX(2*s_elt_a_elt_b_init_loc,
                                           idx_write_init_loc + elt_b_init_loc[elt_b_id]);
          elt_a_elt_b_init_loc = realloc(elt_a_elt_b_init_loc,
                                         sizeof(int) * s_elt_a_elt_b_init_loc * 3);
        }

        for (int j = 0; j < elt_b_init_loc_n[elt_b_id]; j++) {
          elt_a_elt_b_init_loc_stride[elt_a_id]++;
          for (int k = 0; k < 3; k++) {
            elt_a_elt_b_init_loc[3*idx_write_init_loc+k] = elt_b_init_loc[3*elt_b_init_loc_idx[elt_b_id]+k];
          }
          idx_write_init_loc++;
        }

        elt_a_elt_b_n[elt_a_id]++;
        idx_write++;
      }
    }
    // elt_a_elt_b_idx[elt_a_id+1] = idx_write;
    idx_read += n;
  }
  if (idx_write < idx_read) {
      // elt_a_elt_b        = realloc(elt_a_elt_b,        sizeof(int   ) * idx_write);
    elt_a_elt_b_g_num      = realloc(elt_a_elt_b_g_num,      sizeof(PDM_g_num_t) * idx_write);
    // elt_a_elt_b_volume     = realloc(elt_a_elt_b_volume,     sizeof(double     ) * idx_write);
    elt_a_elt_b_init_loc_n = realloc(elt_a_elt_b_init_loc_n, sizeof(int        ) * idx_write);
    elt_a_elt_b_init_loc   = realloc(elt_a_elt_b_init_loc,   sizeof(int        ) * idx_write_init_loc * 3);
  }
  // free(elt_b_init_loc);
  free(elt_b_init_loc_n);
  free(elt_b_init_loc_idx);


  // dbg prints?



  /* Exchange from extracted A to user A */
  PDM_part_to_part_t *ptp_a = NULL;
  PDM_extract_part_part_to_part_get(extrp_mesh_a,
                                    entity_type[0],
                                    &ptp_a,
                                    PDM_OWNERSHIP_KEEP);


  int **user_elt_a_b_n = NULL;
  mi->elt_a_elt_b      = NULL;
  int request_g_num = -1;
  PDM_part_to_part_iexch(ptp_a,
                         comm_kind,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(PDM_g_num_t),
        (const int   **) &elt_a_elt_b_n,
        (const void  **) &elt_a_elt_b_g_num,
                         &user_elt_a_b_n,
        (      void ***) &mi->elt_a_elt_b,
                         &request_g_num);

  PDM_part_to_part_iexch_wait(ptp_a, request_g_num);
  free(elt_a_elt_b_g_num);


  mi->elt_a_elt_b_volume = NULL;
  int request_volume = -1;
  PDM_part_to_part_iexch(ptp_a,
                         comm_kind,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(double),
        (const int   **) &elt_a_elt_b_n,
        (const void  **) &elt_a_elt_b_volume,
                         &user_elt_a_b_n,
        (      void ***) &mi->elt_a_elt_b_volume,
                         &request_volume);

  // exchange faceA_faceB triplets and stride...
  int **user_elt_a_b_init_loc_n = NULL;
  int request_init_loc_n = -1;
  PDM_part_to_part_iexch(ptp_a,
                         comm_kind,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(int),
        (const int   **) &elt_a_elt_b_n,
        (const void  **) &elt_a_elt_b_init_loc_n,
                         &user_elt_a_b_n,
        (      void ***) &user_elt_a_b_init_loc_n,
                         &request_init_loc_n);

  int **user_elt_a_b_init_loc_stride = NULL;
  int **user_elt_a_b_init_loc = NULL;
  int request_init_loc = -1;
  PDM_part_to_part_iexch(ptp_a,
                         comm_kind,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         3*sizeof(int),
        (const int   **) &elt_a_elt_b_init_loc_stride,
        (const void  **) &elt_a_elt_b_init_loc,
                         &user_elt_a_b_init_loc_stride,
        (      void ***) &user_elt_a_b_init_loc,
                         &request_init_loc);

  PDM_part_to_part_iexch_wait(ptp_a, request_volume);
  PDM_part_to_part_iexch_wait(ptp_a, request_init_loc_n);
  PDM_part_to_part_iexch_wait(ptp_a, request_init_loc);
  free(elt_a_elt_b_init_loc_stride);
  free(elt_a_elt_b_init_loc);
  free(elt_a_elt_b_init_loc_n);
  free(elt_a_elt_b_n);


  int  *n_ref_a = NULL;
  int **ref_a   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp_a,
                                 &n_ref_a,
                                 &ref_a);

  int          *user_n_elt_a              = malloc(sizeof(int          ) * mi->n_part_mesh[0]);
  PDM_g_num_t **user_elt_ln_to_gn_a       = malloc(sizeof(PDM_g_num_t *) * mi->n_part_mesh[0]);
  mi->elt_a_elt_b_idx                     = malloc(sizeof(int         *) * mi->n_part_mesh[0]);
  int         **user_elt_a_b_init_loc_idx = malloc(sizeof(int         *) * mi->n_part_mesh[0]);// size = user_a_b_idx[user_n_elt_a]+1
  // int         **user_a_b_init_loc     = malloc(sizeof(int         *) * mi->n_part_mesh[0]);// size = user_a_b_init_loc_idx[user_a_b_idx[user_n_elt_a]] (*3?)
  for (int ipart = 0; ipart < mi->n_part_mesh[0]; ipart++) {
    free(user_elt_a_b_init_loc_stride[ipart]);
    if (mi->mesh_nodal[0] == NULL && mi->mesh[0] != NULL) {
      user_n_elt_a[ipart] = PDM_part_mesh_n_entity_get(mi->mesh[0],
                                                       ipart,
                                                       entity_type[0]);
      PDM_part_mesh_entity_ln_to_gn_get(mi->mesh[0],
                                        ipart,
                                        entity_type[0],
                                        &user_elt_ln_to_gn_a[ipart],
                                        PDM_OWNERSHIP_USER);
    }
    else if (mi->mesh_nodal[0] != NULL) {
      PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_MAX;
      switch (mi->dim_mesh[0]) {
      case 1:
        geom_kind = PDM_GEOMETRY_KIND_RIDGE;
        break;
      case 2:
        geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
        break;
      case 3:
        geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
        break;
      default:
        PDM_error(__FILE__, __LINE__, 0, "invalid dimension %d\n", mi->dim_mesh[0]);
      }

      int n_section = PDM_part_mesh_nodal_n_section_in_geom_kind_get(mi->mesh_nodal[0],
                                                                     geom_kind);
      int *sections_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(mi->mesh_nodal[0],
                                                                          geom_kind);

      user_n_elt_a[ipart] = 0;
      for (int isection = 0; isection < n_section; isection++) {
        int id_section_in_geom_kind = sections_id[isection];
        int id_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(mi->mesh_nodal[0],
                                                                           geom_kind,
                                                                           id_section_in_geom_kind);
        user_n_elt_a[ipart] += PDM_part_mesh_nodal_section_n_elt_get(mi->mesh_nodal[0],
                                                                     id_section,
                                                                     ipart);
      }

      user_elt_ln_to_gn_a[ipart] = malloc(sizeof(PDM_g_num_t) * user_n_elt_a[ipart]);
      for (int isection = 0; isection < n_section; isection++) {
        int id_section_in_geom_kind = sections_id[isection];
        int id_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(mi->mesh_nodal[0],
                                                                           geom_kind,
                                                                           id_section_in_geom_kind);

        int n_elt_section = PDM_part_mesh_nodal_section_n_elt_get(mi->mesh_nodal[0],
                                                                  id_section,
                                                                  ipart);

        PDM_g_num_t *elt_ln_to_gn = PDM_part_mesh_nodal_g_num_get(mi->mesh_nodal[0],
                                                                  id_section,
                                                                  ipart,
                                                                  PDM_OWNERSHIP_KEEP);
        int *parent_num = PDM_part_mesh_nodal_section_parent_num_get(mi->mesh_nodal[0],
                                                                     id_section,
                                                                     ipart,
                                                                     PDM_OWNERSHIP_KEEP);

        for (int ielt = 0; ielt < n_elt_section; ielt++)  {
          int i = ielt;
          if (parent_num != NULL) {
            i = parent_num[ielt];
          }
          user_elt_ln_to_gn_a[ipart][i] = elt_ln_to_gn[ielt];
        }
      }

    }
    else {
      user_n_elt_a[ipart] = 0;
      user_elt_ln_to_gn_a[ipart] = NULL;
    }


    mi->elt_a_elt_b_idx[ipart] = PDM_array_zeros_int(user_n_elt_a[ipart]+1);
    for (int i = 0; i < n_ref_a[ipart]; i++) {
      int elt_a_id = ref_a[ipart][i] - 1;
      mi->elt_a_elt_b_idx[ipart][elt_a_id+1] = user_elt_a_b_n[ipart][i];
    }
    free(user_elt_a_b_n[ipart]);

    for (int i = 0; i < user_n_elt_a[ipart]; i++) {
      mi->elt_a_elt_b_idx[ipart][i+1] += mi->elt_a_elt_b_idx[ipart][i];
    }

    int n = mi->elt_a_elt_b_idx[ipart][user_n_elt_a[ipart]];

    user_elt_a_b_init_loc_idx[ipart] = malloc(sizeof(int) * (n + 1));
    user_elt_a_b_init_loc_idx[ipart][0] = 0;
    int max_init_loc_n = 0;
    for (int i = 0; i < n; i++) {
      max_init_loc_n = PDM_MAX(max_init_loc_n, user_elt_a_b_init_loc_n[ipart][i]);

      user_elt_a_b_init_loc_idx[ipart][i+1] = user_elt_a_b_init_loc_idx[ipart][i] + user_elt_a_b_init_loc_n[ipart][i]*3;
    }

    /* Lexicographic sort on init loc triplets */
    int *order = malloc(sizeof(int) * max_init_loc_n);
    for (int i = 0; i < n; i++) {
      PDM_order_lnum_s(user_elt_a_b_init_loc[ipart] + user_elt_a_b_init_loc_idx[ipart][i],
                       3,
                       order,
                       user_elt_a_b_init_loc_n[ipart][i]);
    }
    free(order);

    free(user_elt_a_b_init_loc_n[ipart]);
  }
  free(user_elt_a_b_n);
  free(user_elt_a_b_init_loc_n);
  free(user_elt_a_b_init_loc_stride);


  int *user_n_elt_b = malloc(sizeof(int) * mi->n_part_mesh[1]);
  for (int ipart = 0; ipart < mi->n_part_mesh[1]; ipart++) {
    if (mi->mesh_nodal[1] == NULL && mi->mesh[1] != NULL) {
      user_n_elt_b[ipart] = PDM_part_mesh_n_entity_get(mi->mesh[1],
                                                       ipart,
                                                       entity_type[1]);
    }
    else if (mi->mesh_nodal[1] != NULL) {
      PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_MAX;
      switch (mi->dim_mesh[1]) {
      case 1:
        geom_kind = PDM_GEOMETRY_KIND_RIDGE;
        break;
      case 2:
        geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
        break;
      case 3:
        geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
        break;
      default:
        PDM_error(__FILE__, __LINE__, 0, "invalid dimension %d\n", mi->dim_mesh[1]);
      }

      int n_section = PDM_part_mesh_nodal_n_section_in_geom_kind_get(mi->mesh_nodal[1],
                                                                     geom_kind);
      int *sections_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(mi->mesh_nodal[1],
                                                                          geom_kind);

      user_n_elt_b[ipart] = 0;
      for (int isection = 0; isection < n_section; isection++) {
        int id_section_in_geom_kind = sections_id[isection];
        int id_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(mi->mesh_nodal[1],
                                                                           geom_kind,
                                                                           id_section_in_geom_kind);
        user_n_elt_b[ipart] += PDM_part_mesh_nodal_section_n_elt_get(mi->mesh_nodal[1],
                                                                     id_section,
                                                                     ipart);
      }
    }
    else {
      user_n_elt_b[ipart] = 0;
    }
  }

  // dbg print?
  mi->ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) user_elt_ln_to_gn_a,
                                                      (const int          *) user_n_elt_a,
                                                                             mi->n_part_mesh[0],
                                                      (const int          *) user_n_elt_b,
                                                                             mi->n_part_mesh[1],
                                                      (const int         **) mi->elt_a_elt_b_idx,       // size = user_n_elt_a+1
                                                      (const int         **) user_elt_a_b_init_loc_idx, // size = user_a_b_idx[user_n_elt_a]+1
                                                      (const int         **) user_elt_a_b_init_loc,     // size = user_a_b_init_loc_idx[user_a_b_idx[user_n_elt_a]] (*3?)
                                                                             mi->comm);
  for (int ipart = 0; ipart < mi->n_part_mesh[0]; ipart++) {
    free(user_elt_a_b_init_loc_idx[ipart]);
    free(user_elt_a_b_init_loc    [ipart]);
    if (mi->mesh_nodal[0] != NULL) {
      free(user_elt_ln_to_gn_a[ipart]);
    }
  }
  free(user_elt_a_b_init_loc_idx);
  free(user_elt_a_b_init_loc    );
  free(user_elt_ln_to_gn_a);

  /* Reverse weights */
  if (1) {
    request_volume = -1;
    PDM_part_to_part_iexch(mi->ptp,
                           comm_kind,
                           PDM_STRIDE_CST_INTERLACED,
                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                           1,
                           sizeof(double),
                           NULL,
          (const void  **) mi->elt_a_elt_b_volume,
                           NULL,
          (void       ***) &mi->elt_b_elt_a_volume,
                           &request_volume);

    PDM_part_to_part_iexch_wait(mi->ptp, request_volume);
  }

  free(user_n_elt_a);
  free(user_n_elt_b);
}










static void
_dump_elementary_vol_vol
(
 PDM_MPI_Comm  comm,
 int           cellA_id,
 int           faceA_id,
 int           cellB_id,
 int           triaA_id,
 int           faceB_id,
 int           triaB_id,
 double       *tetraA_coord,
 double       *triaB_coord
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  char filename[999];
  sprintf(filename, "vol_vol_A%d_%d_%d_B%d_%d_%d_rank_%d.vtk",
          cellA_id,
          faceA_id,
          triaA_id,
          cellB_id,
          faceB_id,
          triaB_id,
          i_rank);


  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 4 + 3);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", tetraA_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", triaB_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  // fprintf(f, "CELLS %d %d\n", 2, 5 + 4);
  // fprintf(f, "4 0 1 2 3\n3 4 5 6\n");

  // fprintf(f, "CELL_TYPES 2\n10\n5\n");

  // fprintf(f, "CELL_DATA 2\n");
  // fprintf(f, "FIELD field 4\n");
  // fprintf(f, "mesh_id 1 2 int\n0 1\n");
  // fprintf(f, "cell_id 1 2 int\n%d %d\n", cellA_id, cellB_id);
  // fprintf(f, "face_id 1 2 int\n%d %d\n", faceA_id, faceB_id);
  // fprintf(f, "tria_id 1 2 int\n%d %d\n", triaA_id, triaB_id);
  fprintf(f, "CELLS %d %d\n", 1 + 3 + 2 + 1, 2 + 3*3 + 2*4 + 4);
  fprintf(f, "1 0\n2 0 1\n2 0 2\n2 0 3\n3 1 2 3\n3 0 1 2\n3 4 5 6\n");

  fprintf(f, "CELL_TYPES 7\n1\n3\n3\n3\n5\n5\n5\n");

  fprintf(f, "CELL_DATA 7\n");
  fprintf(f, "FIELD field 4\n");
  fprintf(f, "mesh_id 1 7 int\n0 0 0 0 0 0 1\n");
  fprintf(f, "cell_id 1 7 int\n%d %d %d %d %d %d %d\n", cellA_id, cellA_id, cellA_id, cellA_id, cellA_id, cellA_id, cellB_id);
  fprintf(f, "face_id 1 7 int\n%d %d %d %d %d %d %d\n", faceA_id, faceA_id, faceA_id, faceA_id, faceA_id, faceA_id, faceB_id);
  fprintf(f, "tria_id 1 7 int\n%d %d %d %d %d %d %d\n", triaA_id, triaA_id, triaA_id, triaA_id, triaA_id, triaA_id, triaB_id);

  fclose(f);
}


static
void
_mesh_intersection_vol_vol
(
 PDM_mesh_intersection_t *mi,
 int                     *a_to_b_idx,
 int                     *a_to_b
)
{
  /* method : 0 -> with pointers, 1 -> without, 2 -> without (more robust) */
  int method = 2;

  int dbg_enabled = 0;


  int          n_vtx        [2] = {0};
  int          n_face       [2] = {0};
  int          n_cell       [2] = {0};
  double      *vtx_coord    [2] = {NULL};
  int         *cell_face_idx[2] = {NULL};
  int         *cell_face    [2] = {NULL};
  int         *face_vtx_idx [2] = {NULL};
  int         *face_vtx     [2] = {NULL};
  PDM_g_num_t *vtx_ln_to_gn [2] = {NULL};
  // PDM_g_num_t *face_ln_to_gn[2] = {NULL};
  PDM_g_num_t *cell_ln_to_gn[2] = {NULL};

  PDM_part_mesh_t *extract_part_mesh[2] = {NULL};

  int          is_owner_face_vtx[2] = {0};

  /* Get extract volume meshes */
  for (int i = 0; i < 2; i++) {
    n_vtx[i] = PDM_extract_part_n_entity_get(mi->extrp_mesh[i],
                                             0,
                                             PDM_MESH_ENTITY_VTX);

    PDM_extract_part_parent_ln_to_gn_get(mi->extrp_mesh[i],
                                         0,
                                         PDM_MESH_ENTITY_VTX,
                                         &vtx_ln_to_gn[i],
                                         PDM_OWNERSHIP_KEEP);

    PDM_extract_part_vtx_coord_get(mi->extrp_mesh[i],
                                   0,
                                   &vtx_coord[i],
                                   PDM_OWNERSHIP_KEEP);

    n_cell[i] = mi->extrp_mesh[i]->n_target[0];
    cell_ln_to_gn[i] = mi->extrp_mesh[i]->target_gnum[0];

    if (dbg_enabled) {
      log_trace("mesh %d is nodal? %d\n", i, mi->mesh[i] == NULL);
      PDM_log_trace_array_long(cell_ln_to_gn[i], n_cell[i], "cell_ln_to_gn : ");
    }

    if (mi->mesh[i] == NULL) {
      /* From part mesh nodal */
      PDM_part_mesh_nodal_elmts_t *extract_pmne = NULL;
      PDM_extract_part_part_mesh_nodal_get(mi->extrp_mesh[i],
                                           &extract_pmne,
                                           PDM_OWNERSHIP_USER);

      PDM_part_mesh_nodal_t *extract_pmn = PDM_part_mesh_nodal_create(mi->dim_mesh[i],
                                                                      1,
                                                                      mi->comm);

      PDM_part_mesh_nodal_coord_set(extract_pmn,
                                    0,
                                    n_vtx       [i],
                                    vtx_coord   [i],
                                    vtx_ln_to_gn[i],
                                    PDM_OWNERSHIP_USER);

      PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(extract_pmn,
                                                    extract_pmne);

      // Convert part_mesh_nodal to part_mesh
      if (dbg_enabled) {
        log_trace(">> PDM_part_mesh_nodal_to_part_mesh\n");

        int i_rank;
        PDM_MPI_Comm_rank(mi->comm, &i_rank);

        // char name[999];
        // sprintf(name, "extract_pmn_vol_mesh%d", i);
        // PDM_part_mesh_nodal_dump_vtk(extract_pmn,
        //                              PDM_GEOMETRY_KIND_VOLUMIC,
        //                              name);

        int n_section    = PDM_part_mesh_nodal_elmts_n_section_get  (extract_pmne);
        int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(extract_pmne);
        for (int isection = 0; isection < n_section; isection++) {

          int id_section = sections_id[isection];
          int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(extract_pmne, id_section, 0);

          PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extract_pmne, id_section);

          char filename[999];
          sprintf(filename, "extrp_pmne_mesh%d_section%d_rank%d.vtk", i, isection, i_rank);

          if (t_elt == PDM_MESH_NODAL_POLY_3D) {

            int          _n_face = 0;
            PDM_g_num_t *_face_ln_to_gn = NULL;
            int         *_face_vtx_idx = NULL;
            int         *_face_vtx = NULL;
            PDM_g_num_t *numabs = NULL;
            int         *_cell_face_idx = NULL;
            int         *_cell_face = NULL;
            int         *parent_num = NULL;
            PDM_g_num_t *parent_entity_g_num = NULL;

            PDM_part_mesh_nodal_elmts_section_poly3d_get(extract_pmne,
                                                         id_section,
                                                         0,
                                                         &_n_face,
                                                         &_face_ln_to_gn,
                                                         &_face_vtx_idx,
                                                         &_face_vtx,
                                                         &numabs,
                                                         &_cell_face_idx,
                                                         &_cell_face,
                                                         &parent_num,
                                                         &parent_entity_g_num,
                                                         PDM_OWNERSHIP_KEEP);
            int *__face_vtx_idx = malloc(sizeof(int) * (_cell_face_idx[n_elt] + 1));
            __face_vtx_idx[0] = 0;
            for (int k = 0; k < _cell_face_idx[n_elt]; k++) {
              int face_id = PDM_ABS(_cell_face[k]) - 1;
              __face_vtx_idx[k+1] = __face_vtx_idx[k] + _face_vtx_idx[face_id+1] - _face_vtx_idx[face_id];
            }

            PDM_g_num_t *__face_cell_ln_to_gn = malloc(sizeof(PDM_g_num_t) * _cell_face_idx[n_elt]);
            int *__face_vtx = malloc(sizeof(int) * __face_vtx_idx[_cell_face_idx[n_elt]]);
            int idx = 0;
            for (int k = 0; k < n_elt; k++) {
              for (int iface = _cell_face_idx[k]; iface < _cell_face_idx[k+1]; iface++) {
                int face_id = PDM_ABS(_cell_face[iface]) - 1;
                __face_cell_ln_to_gn[iface] = parent_entity_g_num[k];
                if (_cell_face[iface] > 0) {
                  for (int j = _face_vtx_idx[face_id]; j < _face_vtx_idx[face_id+1]; j++) {
                    __face_vtx[idx++] = _face_vtx[j];
                  }
                }
                else {
                  for (int j = _face_vtx_idx[face_id+1]-1; j >= _face_vtx_idx[face_id]; j--) {
                    __face_vtx[idx++] = _face_vtx[j];
                  }
                }
              }
            }

            PDM_vtk_write_polydata(filename,
                                   n_vtx[i],
                                   vtx_coord[i],
                                   vtx_ln_to_gn[i],
                                   _cell_face_idx[n_elt],
                                   __face_vtx_idx,
                                   __face_vtx,
                                   __face_cell_ln_to_gn,
                                   NULL);
            free(__face_vtx);
            free(__face_vtx_idx);
            free(__face_cell_ln_to_gn);

          }
          else {

            int         *connec = NULL;
            PDM_g_num_t *numabs = NULL;
            int         *parent_num = NULL;
            PDM_g_num_t *parent_entity_g_num = NULL;
            PDM_part_mesh_nodal_elmts_section_std_get(extract_pmne,
                                                      id_section,
                                                      0,
                                                      &connec,
                                                      &numabs,
                                                      &parent_num,
                                                      &parent_entity_g_num,
                                                      PDM_OWNERSHIP_KEEP);
            PDM_vtk_write_std_elements(filename,
                                       n_vtx[i],
                                       vtx_coord[i],
                                       vtx_ln_to_gn[i],
                                       t_elt,
                                       n_elt,
                                       connec,
                                       parent_entity_g_num,
                                       0, NULL, NULL);
          }
        }
      }
      extract_part_mesh[i] = PDM_part_mesh_nodal_to_part_mesh(extract_pmn,
                                                              PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                                              PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_NONE);

      PDM_part_mesh_nodal_free(extract_pmn);


      PDM_part_mesh_connectivity_get(extract_part_mesh[i],
                                     0,
                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                     &cell_face    [i],
                                     &cell_face_idx[i],
                                     PDM_OWNERSHIP_BAD_VALUE);

      PDM_part_mesh_connectivity_get(extract_part_mesh[i],
                                     0,
                                     PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                     &face_vtx    [i],
                                     &face_vtx_idx[i],
                                     PDM_OWNERSHIP_BAD_VALUE);

      n_face[i] = PDM_part_mesh_n_entity_get(extract_part_mesh[i],
                                             0,
                                             PDM_MESH_ENTITY_FACE);
    }
    else {
      /* From part mesh */
      PDM_extract_part_connectivity_get(mi->extrp_mesh[i],
                                        0,
                                        PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                        &cell_face    [i],
                                        &cell_face_idx[i],
                                        PDM_OWNERSHIP_KEEP);

      n_face[i] = PDM_extract_part_connectivity_get(mi->extrp_mesh[i],
                                                    0,
                                                    PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                    &face_vtx    [i],
                                                    &face_vtx_idx[i],
                                                    PDM_OWNERSHIP_KEEP);

      if (face_vtx[i] == NULL) {
        is_owner_face_vtx[i] = 1;

        int *face_edge = NULL;
        PDM_extract_part_connectivity_get(mi->extrp_mesh[i],
                                          0,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          &face_edge,
                                          &face_vtx_idx[i],
                                          PDM_OWNERSHIP_KEEP);

        int *edge_vtx     = NULL;
        int *edge_vtx_idx = NULL;
        PDM_extract_part_connectivity_get(mi->extrp_mesh[i],
                                          0,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &edge_vtx,
                                          &edge_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);

        PDM_compute_face_vtx_from_face_and_edge(n_face[i],
                                                face_vtx_idx[i],
                                                face_edge,
                                                edge_vtx,
                                                &face_vtx[i]);
      }
    }
  }

  if (dbg_enabled) {
    for (int i = 0; i < n_cell[0]; i++) {
      log_trace(PDM_FMT_G_NUM" (%d faces) -> ",
                cell_ln_to_gn[0][i],
                cell_face_idx[0][i+1] - cell_face_idx[0][i]);
      for (int j = a_to_b_idx[i]; j < a_to_b_idx[i+1]; j++) {
        log_trace(PDM_FMT_G_NUM" ", cell_ln_to_gn[1][a_to_b[j]]);
      }
      log_trace("\n");
    }
  }

  /*
   * Panic vtk
   */
  if (dbg_enabled) {
    for (int i = 0; i < 2; i++) {
      char name[99];
      sprintf(name, "mesh_intersection_vol_vol_mesh%d", i);
      _export_ensight3d(mi->comm,
                        name,
                        n_cell       [i],
                        n_face       [i],
                        n_vtx        [i],
                        cell_face_idx[i],
                        cell_face    [i],
                        face_vtx_idx [i],
                        face_vtx     [i],
                        vtx_coord    [i],
                        cell_ln_to_gn[i],
                        vtx_ln_to_gn [i]);
    }
  }

  /*
   *  Triangulate all faces of mesh B
   */
  int max_face_vtx_n = 0;
  int s_triaA_vtxA = 0;
  for (int i = 0; i < n_face[0]; i++) {
    int face_vtx_n = face_vtx_idx[0][i+1] - face_vtx_idx[0][i];
    max_face_vtx_n = PDM_MAX(max_face_vtx_n, face_vtx_n);
    s_triaA_vtxA   = PDM_MAX(s_triaA_vtxA,   face_vtx_n);
  }
  s_triaA_vtxA = 3*(s_triaA_vtxA - 2);

  int *faceB_triaB_idx = malloc(sizeof(int) * (n_face[1] + 1));
  faceB_triaB_idx[0] = 0;
  for (int i = 0; i < n_face[1]; i++) {
    int face_vtx_n = face_vtx_idx[1][i+1] - face_vtx_idx[1][i];
    max_face_vtx_n = PDM_MAX(max_face_vtx_n, face_vtx_n);
    int face_tria_n = face_vtx_n - 2;
    faceB_triaB_idx[i+1] = faceB_triaB_idx[i] + face_tria_n;
  }

  PDM_triangulate_state_t *tri_state = PDM_triangulate_state_create(max_face_vtx_n);


  int *triaB_vtxB = malloc(sizeof(int) * faceB_triaB_idx[n_face[1]] * 3);
  int *triaA_vtxA = malloc(sizeof(int) * s_triaA_vtxA);

  for (int faceB_id = 0; faceB_id < n_face[1]; faceB_id++) {

    int *_face_vtx = face_vtx[1] + face_vtx_idx[1][faceB_id];
    int face_vtx_n = face_vtx_idx[1][faceB_id+1] - face_vtx_idx[1][faceB_id];

    int *_tria_vtx  = triaB_vtxB + 3*faceB_triaB_idx[faceB_id];

    int n_tria;
    if (face_vtx_n == 3) {
      /* Triangular face */
      n_tria = 1;
      memcpy(_tria_vtx, _face_vtx, sizeof(int) * 3);
    }
    else if (face_vtx_n == 4) {
      /* Quadrilateral face */
      n_tria = PDM_triangulate_quadrangle(3,
                                          vtx_coord[1],
                                          NULL,
                                          _face_vtx,
                                          _tria_vtx);
    }
    else {
      /* Polygonal face */
      n_tria = PDM_triangulate_polygon(3,
                                       face_vtx_n,
                                       vtx_coord[1],
                                       NULL,
                                       _face_vtx,
                                       PDM_TRIANGULATE_MESH_DEF,
                                       _tria_vtx,
                                       tri_state);
    }

    assert(n_tria == faceB_triaB_idx[faceB_id+1] - faceB_triaB_idx[faceB_id]);

  } // End of loop on faces B



  double *cellA_volume = malloc(sizeof(double) * n_cell[0]);
  double *cellA_center = malloc(sizeof(double) * n_cell[0] * 3);
  PDM_geom_elem_polyhedra_properties_triangulated(1,
                                                  n_cell       [0],
                                                  n_face       [0],
                                                  face_vtx_idx [0],
                                                  face_vtx     [0],
                                                  cell_face_idx[0],
                                                  cell_face    [0],
                                                  n_vtx        [0],
                                                  vtx_coord    [0],
                                                  cellA_volume,
                                                  cellA_center,
                                                  NULL,
                                                  NULL);

  if (1) {
    for (int i = 0; i < n_cell[0]; i++) {
      if (cellA_volume[i] < 0) {
        printf("!! cell A "PDM_FMT_G_NUM" : volume = %e\n", cell_ln_to_gn[0][i], cellA_volume[i]);
        // log_trace("flip cellA %d (%d faces)\n", i, cell_face_idx[0][i+1] - cell_face_idx[0][i]);
        // for (int j = cell_face_idx[0][i]; j < cell_face_idx[0][i+1]; j++) {
        //   cell_face[0][j] = -cell_face[0][j];
        // }
      }
    }
  }

  double *a_to_b_volume = malloc(sizeof(double) * a_to_b_idx[n_cell[0]]);

  double tetraA_coord[12];
  double triaB_coord[9];

  // dbg_enabled = 0;
  for (int cellA_id = 0; cellA_id < n_cell[0]; cellA_id++) {

    /* Compute a 'center' point to tetrahedrize current cell A */
    // Reference vertex = 1st vertex of first face
    // Use cell center instead????
    int ref_vtxA_id = -1;
    int first_face_id = PDM_ABS(cell_face[0][cell_face_idx[0][cellA_id]]) - 1;

    if (mi->tetraisation_pt_type == 0) {
      ref_vtxA_id = face_vtx[0][face_vtx_idx[0][first_face_id]];
    }


    if (dbg_enabled) {
      log_trace("cellA_id %d ("PDM_FMT_G_NUM") : ref_vtxA_id = %d\n",
                cellA_id, cell_ln_to_gn[0][cellA_id], ref_vtxA_id);
    }

    if (mi->tetraisation_pt_type == 0) {
      memcpy(&tetraA_coord[0], &vtx_coord[0][3*(ref_vtxA_id-1)], sizeof(double)*3);
    }

    else if (mi->tetraisation_pt_type == 1) {
      memcpy(&tetraA_coord[0], &cellA_center[3*cellA_id], sizeof(double)*3);
    }

    else if (mi->tetraisation_pt_type == 2) {
      tetraA_coord[0] = mi->tetraisation_pt_coord[0];
      tetraA_coord[1] = mi->tetraisation_pt_coord[1];
      tetraA_coord[2] = mi->tetraisation_pt_coord[2];
    }

    /* Initialize intersection volumes to zero */
    for (int icellB = a_to_b_idx[cellA_id]; icellB < a_to_b_idx[cellA_id+1]; icellB++) {
      a_to_b_volume[icellB] = 0.;
    }

    /* Loop on faces A */
    for (int ifaceA = cell_face_idx[0][cellA_id]; ifaceA < cell_face_idx[0][cellA_id+1]; ifaceA++) {

      int faceA_id   = PDM_ABS (cell_face[0][ifaceA]) - 1;
      int faceA_sign = PDM_SIGN(cell_face[0][ifaceA]);

       if (dbg_enabled) {
        log_trace("  faceA_id %d, sign = %d\n", faceA_id, faceA_sign);
      }

      /* Triangulate current face A */
      int *_face_vtx = face_vtx[0] + face_vtx_idx[0][faceA_id];
      int face_vtx_n = face_vtx_idx[0][faceA_id+1] - face_vtx_idx[0][faceA_id];

      int n_tria;
      if (face_vtx_n == 3) {
        /* Triangular face */
        n_tria = 1;
        memcpy(triaA_vtxA, _face_vtx, sizeof(int) * 3);
      }
      else if (face_vtx_n == 4) {
        /* Quadrilateral face */
        n_tria = PDM_triangulate_quadrangle(3,
                                            vtx_coord[0],
                                            NULL,
                                            _face_vtx,
                                            triaA_vtxA);
      }
      else {
        /* Polygonal face */
        n_tria = PDM_triangulate_polygon(3,
                                         face_vtx_n,
                                         vtx_coord[0],
                                         NULL,
                                         _face_vtx,
                                         PDM_TRIANGULATE_MESH_DEF,
                                         triaA_vtxA,
                                         tri_state);
      }

      /* Loop on triangles A */
      for (int triaA_id = 0; triaA_id < n_tria; triaA_id++) {

        int *_triaA_vtxA = triaA_vtxA + 3*triaA_id;
        if (dbg_enabled) {
          PDM_log_trace_array_int(_triaA_vtxA, 3, "    _triaA_vtxA : ");
        }

        /* Ignore if current triangle contains reference vertex */
        if (_triaA_vtxA[0] == ref_vtxA_id ||
            _triaA_vtxA[1] == ref_vtxA_id ||
            _triaA_vtxA[2] == ref_vtxA_id) {
            continue;
        }

        if (faceA_sign > 0) {
          for (int ivtx = 0; ivtx < 3; ivtx++) {
            int vtxA_id = _triaA_vtxA[ivtx] - 1;
            memcpy(&tetraA_coord[3*(ivtx+1)], &vtx_coord[0][3*vtxA_id], sizeof(double)*3);
          }
        }
        else {
          for (int ivtx = 0; ivtx < 3; ivtx++) {
            int vtxA_id = _triaA_vtxA[2-ivtx] - 1;
            memcpy(&tetraA_coord[3*(ivtx+1)], &vtx_coord[0][3*vtxA_id], sizeof(double)*3);
          }
        }

        double mat[3][3];
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            mat[i][j] = tetraA_coord[3*(j+1)+i] - tetraA_coord[i];
          }
        }

        double det = mat[0][0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])
        -            mat[1][0]*(mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2])
        +            mat[2][0]*(mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2]);

        // log_trace("det = %f\n", det);

        PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
        if (det == 0.) {
          continue;
        }
        PDM_GCC_SUPPRESS_WARNING_POP

        double idet = 1./det;


        /* Loop on candidate cells B */
        for (int icellB = a_to_b_idx[cellA_id]; icellB < a_to_b_idx[cellA_id+1]; icellB++) {

          int cellB_id = a_to_b[icellB];

          if (dbg_enabled) {
            log_trace("      icellB %d, cellB_id %d ("PDM_FMT_G_NUM")\n",
                      icellB, cellB_id, cell_ln_to_gn[1][cellB_id]);
          }

          /* Loop on faces B */
          for (int ifaceB = cell_face_idx[1][cellB_id]; ifaceB < cell_face_idx[1][cellB_id+1]; ifaceB++) {

            int faceB_id   = PDM_ABS (cell_face[1][ifaceB]) - 1;
            int faceB_sign = PDM_SIGN(cell_face[1][ifaceB]);

            if (dbg_enabled) {
              log_trace("        faceB_id %d, sign = %d\n", faceB_id, faceB_sign);
            }

            /* Loop on triangles B */
            for (int triaB_id = faceB_triaB_idx[faceB_id]; triaB_id < faceB_triaB_idx[faceB_id+1]; triaB_id++) {

              int *_triaB_vtxB = triaB_vtxB + 3*triaB_id;

              if (dbg_enabled) {
                log_trace("          triaB_id %d : ", triaB_id);
                PDM_log_trace_array_int(_triaB_vtxB, 3, "");
              }

              if (faceB_sign > 0) {
                for (int ivtx = 0; ivtx < 3; ivtx++) {
                  int vtxB_id = _triaB_vtxB[ivtx] - 1;
                  memcpy(&triaB_coord[3*ivtx], &vtx_coord[1][3*vtxB_id], sizeof(double)*3);
                }
              }
              else {
                for (int ivtx = 0; ivtx < 3; ivtx++) {
                  int vtxB_id = _triaB_vtxB[2-ivtx] - 1;
                  memcpy(&triaB_coord[3*ivtx], &vtx_coord[1][3*vtxB_id], sizeof(double)*3);
                }
              }

              if (dbg_enabled && 0) {
                _dump_elementary_vol_vol(mi->comm,
                                         cellA_id,
                                         faceA_id,
                                         cellB_id,
                                         triaA_id,
                                         faceB_id,
                                         triaB_id,
                                         tetraA_coord,
                                         triaB_coord);
              }


              /* Transform triangle to current tetra's local frame */
              for (int ivtx = 0; ivtx < 3; ivtx++) {
                double rhs[3];
                for (int i = 0; i < 3; i++) {
                  rhs[i] = triaB_coord[3*ivtx + i] - tetraA_coord[i];
                }

                triaB_coord[3*ivtx    ] = idet *
                (  rhs[0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])
                 - rhs[1]*(mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2])
                 + rhs[2]*(mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2]));

                triaB_coord[3*ivtx + 1] = idet *
                (  mat[0][0]*(rhs[1]*mat[2][2] - rhs[2]*mat[1][2])
                 - mat[1][0]*(rhs[0]*mat[2][2] - rhs[2]*mat[0][2])
                 + mat[2][0]*(rhs[0]*mat[1][2] - rhs[1]*mat[0][2]));

                triaB_coord[3*ivtx + 2] = idet *
                (  mat[0][0]*(mat[1][1]*rhs[2] - mat[2][1]*rhs[1])
                 - mat[1][0]*(mat[0][1]*rhs[2] - mat[2][1]*rhs[0])
                 + mat[2][0]*(mat[0][1]*rhs[1] - mat[1][1]*rhs[0]));
              }

              if (dbg_enabled && 1) {
                // check inverse transform
                for (int i = 0; i < 3; i++) {
                  int vtxB_id = _triaB_vtxB[i] - 1;
                  if (faceB_sign < 0) {
                    vtxB_id = _triaB_vtxB[2-i] - 1;
                  }
                  double err = 0;
                  for (int j = 0; j < 3; j++) {
                    double u = triaB_coord[3*i  ];
                    double v = triaB_coord[3*i+1];
                    double w = triaB_coord[3*i+2];
                    double delta = vtx_coord[1][3*vtxB_id+j] -
                    ((1-u-v-w) * tetraA_coord[  j] +
                     u         * tetraA_coord[3+j] +
                     v         * tetraA_coord[6+j] +
                     w         * tetraA_coord[9+j]);
                    err += delta*delta;
                  }

                  if (err > 1e-12) {
                    log_trace("!!! error = %e\n", sqrt(err));
                  }
                }
              }


              /* Perform elementary computation */

              double volume = 0;

              if (method == 0) {
                int     vtk_n_faceA      = 0;
                int     vtk_n_faceB      = 0;
                double *local_vtx_coordA = NULL;
                int     local_n_vtxA     = 0;
                int    *local_face_vtxA  = NULL;
                double *local_vtx_coordB = NULL;
                int     local_n_vtxB     = 0;
                int    *local_face_vtxB  = NULL;
                volume = PDM_mesh_intersection_vol_vol_atomic_compute(triaB_coord,
                                                                      &local_vtx_coordA,
                                                                      &local_n_vtxA,
                                                                      &local_face_vtxA,
                                                                      &vtk_n_faceA,
                                                                      &local_vtx_coordB,
                                                                      &local_n_vtxB,
                                                                      &local_face_vtxB,
                                                                      &vtk_n_faceB);
                free(local_vtx_coordA);
                free(local_face_vtxA);
                free(local_vtx_coordB);
                free(local_face_vtxB);
              }
              else if (method == 1) {
                volume = PDM_mesh_intersection_vol_vol_atomic_compute2(triaB_coord);
              }
              else {
                volume = PDM_mesh_intersection_vol_vol_atomic_compute3(triaB_coord);
              }

              if (dbg_enabled) {
                log_trace("            volume elem = %20.16f\n", volume);
                if (1) {//volume != 0) {
                  log_trace("********\n");
                  for (int i = 0; i < 3; i++) {
                    log_trace("double pt%d[3] = {%21.17e,%21.17e,%21.17e};\n",
                              i, triaB_coord[3*i], triaB_coord[3*i+1], triaB_coord[3*i+2]);
                  }
                  log_trace("********\n");
                }
              }

              /* Add elementray volume contribution */
              // log_trace("sign(det) = %d / faceA_sign = %d\n", PDM_SIGN(det), faceA_sign);
              a_to_b_volume[icellB] += volume * PDM_ABS(det);

            } // End of loop on triangles of current face B

          } // End of loop on faces of current cell B

        } // End of loop on candidate cells B for current cell A

      } // End of loop on triangles of current face A

    } // End of loop on faces of current cell A

    if (dbg_enabled) {
      for (int icellB = a_to_b_idx[cellA_id]; icellB < a_to_b_idx[cellA_id+1]; icellB++) {
        int cellB_id = a_to_b[icellB];
        if (a_to_b_volume[icellB] > 0.00001) {log_trace("---> BIG");}
        log_trace("cellA %d ("PDM_FMT_G_NUM") cellB %d ("PDM_FMT_G_NUM"), volume = %20.16f\n",
                  cellA_id, cell_ln_to_gn[0][cellA_id],
                  cellB_id, cell_ln_to_gn[1][cellB_id],
                  a_to_b_volume[icellB]);
      }
    }

  } // End of loop on cells A

  free(faceB_triaB_idx);
  free(triaB_vtxB);
  free(triaA_vtxA);
  PDM_triangulate_state_destroy(tri_state);

  if (dbg_enabled) {
    // Crude check
    double l_total_volume_AB = 0;
    for (int i = 0; i < a_to_b_idx[n_cell[0]]; i++) {
      l_total_volume_AB += a_to_b_volume[i];
    }

    double l_total_volume_A  = 0;
    // double *cellA_volume = malloc(sizeof(double) * n_cellA);
    // double *cellA_center = malloc(sizeof(double) * n_cellA * 3);
    // PDM_geom_elem_polyhedra_properties_triangulated(1,
    //                                                 n_cell       [0],
    //                                                 n_face       [0],
    //                                                 face_vtx_idx [0],
    //                                                 face_vtx     [0],
    //                                                 cell_face_idx[0],
    //                                                 cell_face    [0],
    //                                                 n_vtx        [0],
    //                                                 vtx_coord    [0],
    //                                                 cellA_volume,
    //                                                 cellA_center,
    //                                                 NULL,
    //                                                 NULL);
    // free(cellA_center);

    for (int cellA_id = 0; cellA_id < n_cell[0]; cellA_id++) {
      l_total_volume_A += cellA_volume[cellA_id];
    }


    double g_total_volume_AB;
    PDM_MPI_Allreduce(&l_total_volume_AB, &g_total_volume_AB, 1,
                      PDM_MPI_DOUBLE, PDM_MPI_SUM, mi->comm);

    double g_total_volume_A;
    PDM_MPI_Allreduce(&l_total_volume_A, &g_total_volume_A, 1,
                      PDM_MPI_DOUBLE, PDM_MPI_SUM, mi->comm);

    log_trace("total volume of A inter B : local = %20.16f, global = %20.16f (%3.3f%%)\n",
              l_total_volume_AB, g_total_volume_AB,
              100*g_total_volume_AB / g_total_volume_A);

    // cas cube, translation (0.5,0.5,0.5)
    // double exact = 0.5;
    // log_trace("error : absolute = %e, relative = %e\n",
    //           PDM_ABS(g_total_volume_AB - exact),
    //           PDM_ABS(g_total_volume_AB - exact)/exact);

    // debug
    mi->local_vol_A_B  = l_total_volume_AB;
    mi->global_vol_A_B = g_total_volume_AB;
    mi->global_vol_A   = g_total_volume_A;
  }

  free(cellA_center);
  free(cellA_volume);

  for (int i = 0; i < 2; i++) {
    if (is_owner_face_vtx[i]) {
      free(face_vtx[i]);
    }

    if (extract_part_mesh[i] != NULL) {
      PDM_part_mesh_free(extract_part_mesh[i]);
    }
  }

  _build_ptp(mi,
             a_to_b_idx,
             a_to_b,
             a_to_b_volume);

  free(a_to_b_volume);
}

static
void
_mesh_intersection_vol_surf
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  PDM_UNUSED(mi);
  PDM_UNUSED(extrp_mesh_a);
  PDM_UNUSED(extrp_mesh_b);
  PDM_UNUSED(redistribute_box_a_to_box_b_idx);
  PDM_UNUSED(redistribute_box_a_to_box_b);

  /*
   * Panic vtk
   */
  if(0 == 1) {
    _export_vtk_3d("extrp_mesh_a", extrp_mesh_a);
    _export_vtk_2d("extrp_mesh_b", extrp_mesh_b);
  }

}

static PDM_polygon_status_t
_intersect_ray_face
(
       double                  *face_center,
       double                  *face_normal,
       double                  *face_coord,
       double                  *vtx_coord,
 const int                      face_vtx_n,
 const int                     *face_vtx,
 const double                  *ray_origin,
 const double                  *ray_direction,
       double                  *intersection_coord
 )
{
  double face_bound[6] = {
    HUGE_VAL, -HUGE_VAL,
    HUGE_VAL, -HUGE_VAL,
    HUGE_VAL, -HUGE_VAL
  };
  for (int i = 0; i < face_vtx_n; i++) {
    int vtx_id = face_vtx[i] - 1;
    double *vc = vtx_coord + 3*vtx_id;
    for (int j = 0; j < 3; j++) {
      face_coord[3*i+j] = vc[j];
      face_bound[2*j  ] = PDM_MIN(face_bound[2*j  ], vc[j]);
      face_bound[2*j+1] = PDM_MAX(face_bound[2*j+1], vc[j]);
    }
  }

  // Inflate the face's bounding box
  double d = 0.;
  for (int j = 0; j < 3; j++) {
    d += (face_bound[2*j+1] - face_bound[2*j])*(face_bound[2*j+1] - face_bound[2*j]);
  }
  d = 0.1*sqrt(d);
  for (int j = 0; j < 3; j++) {
    face_bound[2*j  ] -= d;
    face_bound[2*j+1] += d;
  }

  double t;
  PDM_polygon_status_t stat = PDM_polygon_ray_intersection(ray_origin,
                                                           ray_direction,
                                                           face_vtx_n,
                                                           face_coord,
                                                           face_center,
                                                           face_normal,
                                                           face_bound,
                                                           intersection_coord,
                                                           &t,
                                                           NULL);
  if (stat == PDM_POLYGON_INSIDE) {
    if (t < 0 || t > 1) {
      stat = PDM_POLYGON_OUTSIDE;
    }
  }

  return stat;
}

static
void
_mesh_intersection_vol_line
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  int i_rank;
  PDM_MPI_Comm_rank(mi->comm, &i_rank);

  int dbg_enabled = 0;

  int *cellA_lineB_idx = redistribute_box_a_to_box_b_idx;
  int *cellA_lineB     = redistribute_box_a_to_box_b;

  /* Get connectivities and coordinates */
  int          n_cellA = 0;
  int          n_faceA = 0;
  int          n_vtxA  = 0;
  int         *cellA_faceA_idx = NULL;
  int         *cellA_faceA     = NULL;
  int         *faceA_vtxA_idx  = NULL;
  int         *faceA_vtxA      = NULL;
  double      *vtxA_coord      = NULL;
  PDM_g_num_t *cellA_ln_to_gn  = NULL;
  PDM_g_num_t *vtxA_ln_to_gn   = NULL;
  PDM_ownership_t owner_face_vtxA = _get_extracted_mesh_vol(extrp_mesh_a,
                                                            &n_cellA,
                                                            &n_faceA,
                                                            &n_vtxA,
                                                            &cellA_faceA_idx,
                                                            &cellA_faceA,
                                                            &faceA_vtxA_idx,
                                                            &faceA_vtxA,
                                                            &vtxA_coord,
                                                            &cellA_ln_to_gn,
                                                            &vtxA_ln_to_gn);

  double *vtx_coordB = NULL;
  PDM_extract_part_vtx_coord_get(extrp_mesh_b, 0, &vtx_coordB, PDM_OWNERSHIP_KEEP);

  int  *edgeB_vtxB      = NULL;
  int  *edgeB_vtxB_idx  = NULL;
  PDM_extract_part_connectivity_get(extrp_mesh_b, 0, PDM_CONNECTIVITY_TYPE_EDGE_VTX , &edgeB_vtxB , &edgeB_vtxB_idx , PDM_OWNERSHIP_KEEP);


  if (dbg_enabled) {
    // _export_vtk_3d("extrp_mesh_a", extrp_mesh_a);
    // _export_vtk_3d("extrp_mesh_b", extrp_mesh_b);

    _export_ensight3d(mi->comm,
                      "vol_vol_meshA",
                      n_cellA,
                      n_faceA,
                      n_vtxA,
                      cellA_faceA_idx,
                      cellA_faceA,
                      faceA_vtxA_idx,
                      faceA_vtxA,
                      vtxA_coord,
                      cellA_ln_to_gn,
                      vtxA_ln_to_gn);
    _export_vtk_1d("extrp_mesh_b", extrp_mesh_b);

    PDM_log_trace_connectivity_int(cellA_lineB_idx,
                                   cellA_lineB,
                                   n_cellA,
                                   "cellA_lineB : ");
  }

  /*
   * Compute face normal onces
   */
  double *face_normal = malloc(3 * n_faceA * sizeof(double));
  double *face_center = malloc(3 * n_faceA * sizeof(double));

  PDM_geom_elem_polygon_properties(n_faceA,
                                   faceA_vtxA_idx,
                                   faceA_vtxA,
                                   vtxA_coord,
                                   face_normal,
                                   face_center,
                                   NULL,
                                   NULL);

  /*
   * Get max size of vtx per face
   */
  int n_face_vtx_max = 0;
  for(int i_face = 0; i_face < n_faceA; ++i_face) {
    n_face_vtx_max = PDM_MAX(n_face_vtx_max, faceA_vtxA_idx[i_face+1] - faceA_vtxA_idx[i_face]);
  }
  int n_cell_face_max = 0;
  for(int i_cell = 0; i_cell < n_cellA; ++i_cell) {
    n_cell_face_max = PDM_MAX(n_cell_face_max, cellA_faceA_idx[i_cell+1] - cellA_faceA_idx[i_cell]);
  }
  double *poly_coord         = malloc(3 * n_face_vtx_max  * sizeof(double));
  double *intersection_coord = malloc(3 * n_cell_face_max * sizeof(double));
  int    *intersection_stat  = malloc(    n_cell_face_max * sizeof(int   ));

  int         *cellA_lineB_post_idx = malloc((n_cellA+1)              * sizeof(int));
  int         *cellA_lineB_post_n   = malloc((n_cellA)                * sizeof(int));
  int         *cellA_lineB_post     = malloc(cellA_lineB_idx[n_cellA] * sizeof(int));

  /*
   * For each cells we sseek intersection of lines with one faces
   */
  cellA_lineB_post_idx[0] = 0;
  for(int i_cell = 0; i_cell < n_cellA; ++i_cell) {
    cellA_lineB_post_idx[i_cell+1] = cellA_lineB_post_idx[i_cell];
    cellA_lineB_post_n  [i_cell] = 0;
    for(int idx_line = cellA_lineB_idx[i_cell]; idx_line < cellA_lineB_idx[i_cell+1]; ++idx_line) {
      int i_line = cellA_lineB[idx_line];

      int i_vtx1 = edgeB_vtxB[2*i_line  ]-1;
      int i_vtx2 = edgeB_vtxB[2*i_line+1]-1;

      double ray_direction[3] = {
        vtx_coordB[3*i_vtx2  ] - vtx_coordB[3*i_vtx1  ],
        vtx_coordB[3*i_vtx2+1] - vtx_coordB[3*i_vtx1+1],
        vtx_coordB[3*i_vtx2+2] - vtx_coordB[3*i_vtx1+2],
      };

      double ray_origin[3] = {
        vtx_coordB[3*i_vtx1  ],
        vtx_coordB[3*i_vtx1+1],
        vtx_coordB[3*i_vtx1+2],
      };

      int n_intersect = 0;
      int lface = 0;
      for(int idx_face = cellA_faceA_idx[i_cell]; idx_face < cellA_faceA_idx[i_cell+1]; idx_face++) {
        int i_face = PDM_ABS(cellA_faceA[idx_face])-1;

        int *_face_vtx = faceA_vtxA + faceA_vtxA_idx[i_face];
        int face_vtx_n = faceA_vtxA_idx[i_face+1] - faceA_vtxA_idx[i_face];
        intersection_stat[lface] = _intersect_ray_face(&face_center[3*i_face],
                                                       &face_normal[3*i_face],
                                                       poly_coord,
                                                       vtxA_coord,
                                                       face_vtx_n,
                                                       _face_vtx,
                                                       ray_origin,
                                                       ray_direction,
                                                       &intersection_coord[3*lface]);
        // printf("intersection_stat[%i] = %i \n", lface, intersection_stat[lface]);
        if (intersection_stat[lface] == PDM_POLYGON_INSIDE) {
          n_intersect++;
        }
        lface++;
      } /* End face_vtx loop */

      /* Check if ray is purely inside the cell -> localisation of vtx ? */
      // printf("i_cell = %i | i_line = %i | n_intersect = %i \n", i_cell, i_line, n_intersect);

      /* Post-treatment */
      lface = 0;
      for(int idx_face = cellA_faceA_idx[i_cell]; idx_face < cellA_faceA_idx[i_cell+1]; idx_face++) {
        if(intersection_stat[lface] == PDM_POLYGON_INSIDE) {
          cellA_lineB_post_n[i_cell]++;
          cellA_lineB_post[cellA_lineB_post_idx[i_cell+1]++] = i_line;
          break;
        }
        lface++;
      }

    }
  }

  if(0 == 1) {
    // PDM_log_trace_array_long(cellA_ln_to_gn, n_cellA, "cellA_ln_to_gn ::");
    PDM_log_trace_connectivity_int(cellA_lineB_post_idx,
                                   cellA_lineB_post,
                                   n_cellA,
                                   "cellA_lineB_post : ");
  }


  if (owner_face_vtxA == PDM_OWNERSHIP_USER) {
    free(faceA_vtxA);
  }
  free(poly_coord);
  free(face_normal);
  free(face_center);
  free(intersection_coord);
  free(intersection_stat);

  /*
   * Creation du part_to_part
   * may not work if multiple init locations...
   */
  PDM_g_num_t *elt_b_ln_to_gn = extrp_mesh_b->target_gnum[0];
  int         *elt_b_init_loc = extrp_mesh_b->target_location[0];

  int         *elt_a_elt_b_init_loc   = malloc(3 * cellA_lineB_post_idx[n_cellA] * sizeof(int        ));
  PDM_g_num_t *cellA_lineB_post_g_num = malloc(    cellA_lineB_post_idx[n_cellA] * sizeof(PDM_g_num_t));
  // int         *elt_a_elt_b_init_loc_n = malloc(    cellA_lineB_post_idx[n_cellA] * sizeof(int        ));

  // int n_init_loc_tot = 0;
  for(int i_cell = 0; i_cell < n_cellA; ++i_cell) {
    for(int idx_line = cellA_lineB_post_idx[i_cell]; idx_line < cellA_lineB_post_idx[i_cell+1]; ++idx_line) {
      int i_line = cellA_lineB_post[idx_line];
      cellA_lineB_post_g_num[  idx_line  ] = elt_b_ln_to_gn[  i_line  ];
      elt_a_elt_b_init_loc  [3*idx_line  ] = elt_b_init_loc[3*i_line  ];
      elt_a_elt_b_init_loc  [3*idx_line+1] = elt_b_init_loc[3*i_line+1];
      elt_a_elt_b_init_loc  [3*idx_line+2] = elt_b_init_loc[3*i_line+2];

    }
  }

  // int         *elt_a_elt_b_init_loc_n = malloc(    cellA_lineB_post_idx[n_cellA] * sizeof(int        ));
  // for(int i_cell = 0; i_cell < n_cellA; ++i_cell) {
  //   for(int idx_line = cellA_lineB_post_idx[i_cell]; idx_line < cellA_lineB_post_idx[i_cell+1]; ++idx_line) {

  //   }
  // }

  /* Exchange from extracted A to user A */
  PDM_part_to_part_t *ptp_a = NULL;
  PDM_extract_part_part_to_part_get(extrp_mesh_a,
                                    PDM_MESH_ENTITY_CELL,
                                    &ptp_a,
                                    PDM_OWNERSHIP_KEEP);


  PDM_mpi_comm_kind_t comm_kind = PDM_MPI_COMM_KIND_P2P;//COLLECTIVE;
  int **user_elt_a_b_n = NULL;
  mi->elt_a_elt_b      = NULL;
  int request_g_num = -1;
  PDM_part_to_part_iexch(ptp_a,
                         comm_kind,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(PDM_g_num_t),
        (const int   **) &cellA_lineB_post_n,
        (const void  **) &cellA_lineB_post_g_num,
                         &user_elt_a_b_n,
        (      void ***) &mi->elt_a_elt_b,
                         &request_g_num);
  PDM_part_to_part_iexch_wait(ptp_a, request_g_num);

  int **user_elt_a_b_init_loc = NULL;
  int request_init_loc = -1;
  PDM_part_to_part_iexch(ptp_a,
                         comm_kind,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         3 * sizeof(int),
        (const int   **) &cellA_lineB_post_n,
        (const void  **) &elt_a_elt_b_init_loc,
                         &user_elt_a_b_n,
        (      void ***) &user_elt_a_b_init_loc,
                         &request_init_loc);
  PDM_part_to_part_iexch_wait(ptp_a, request_init_loc);

  free(elt_a_elt_b_init_loc);
  free(cellA_lineB_post_g_num);

  int  *n_ref_a = NULL;
  int **ref_a   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp_a,
                                 &n_ref_a,
                                 &ref_a);

  int          *user_n_elt_a              = malloc(mi->n_part_mesh[0] * sizeof(int          ));
  int         **user_elt_a_b_init_loc_idx = malloc(mi->n_part_mesh[0] * sizeof(int         *));
  PDM_g_num_t **user_elt_ln_to_gn_a       = malloc(mi->n_part_mesh[0] * sizeof(PDM_g_num_t *));
  mi->elt_a_elt_b_idx                     = malloc(mi->n_part_mesh[0] * sizeof(int         *));

  /*
   * TODO : merge with init_location
   */

  for (int i_part = 0; i_part < mi->n_part_mesh[0]; i_part++) {
    user_n_elt_a[i_part] = PDM_part_mesh_n_entity_get(mi->mesh[0],
                                                      i_part,
                                                      PDM_MESH_ENTITY_CELL);

    PDM_part_mesh_entity_ln_to_gn_get(mi->mesh[0],
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &user_elt_ln_to_gn_a[i_part],
                                      PDM_OWNERSHIP_USER);

    mi->elt_a_elt_b_idx[i_part] = PDM_array_zeros_int(user_n_elt_a[i_part]+1);

    for (int i = 0; i < n_ref_a[i_part]; i++) {
      int elt_a_id = ref_a[i_part][i] - 1;
      mi->elt_a_elt_b_idx[i_part][elt_a_id+1] = user_elt_a_b_n[i_part][i];
    }
    free(user_elt_a_b_n[i_part]);

    for (int i = 0; i < user_n_elt_a[i_part]; i++) {
      mi->elt_a_elt_b_idx[i_part][i+1] += mi->elt_a_elt_b_idx[i_part][i];
    }

    int n_elt_a_elt_b = mi->elt_a_elt_b_idx[i_part][user_n_elt_a[i_part]];
    user_elt_a_b_init_loc_idx[i_part] = malloc((n_elt_a_elt_b+1) * sizeof(int));
    for(int i = 0; i < n_elt_a_elt_b+1; ++i) {
      user_elt_a_b_init_loc_idx[i_part][i] = 3*i;
    }

    if(0 == 1) {
      PDM_log_trace_connectivity_long(mi->elt_a_elt_b_idx[i_part], mi->elt_a_elt_b[i_part], user_n_elt_a[i_part], "elt_a_elt_b ::");
    }

    // free(user_elt_a_b_idx);
  }
  free(user_elt_a_b_n);

  int *user_n_elt_b = malloc(mi->n_part_mesh[1] * sizeof(int));
  for (int i_part = 0; i_part < mi->n_part_mesh[1]; i_part++) {
    if (mi->mesh_nodal[1] == NULL && mi->mesh[1] != NULL) {
      user_n_elt_b[i_part] = PDM_part_mesh_n_entity_get(mi->mesh[1],
                                                       i_part,
                                                       PDM_MESH_ENTITY_EDGE);
    } else {
      PDM_error(__FILE__, __LINE__, 0, "invalid mesh type \n");
      abort();
    }

  }


  // dbg print?
  mi->ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) user_elt_ln_to_gn_a,
                                                      (const int          *) user_n_elt_a,
                                                                             mi->n_part_mesh[0],
                                                      (const int          *) user_n_elt_b,
                                                                             mi->n_part_mesh[1],
                                                      (const int         **) mi->elt_a_elt_b_idx,       // size = user_n_elt_a+1
                                                      (const int         **) user_elt_a_b_init_loc_idx, // size = user_a_b_idx[user_n_elt_a]+1
                                                      (const int         **) user_elt_a_b_init_loc,     // size = user_a_b_init_loc_idx[user_a_b_idx[user_n_elt_a]] (*3?)
                                                                             mi->comm);

  for (int ipart = 0; ipart < mi->n_part_mesh[0]; ipart++) {
    free(user_elt_a_b_init_loc_idx[ipart]);
    free(user_elt_a_b_init_loc    [ipart]);
    if (mi->mesh_nodal[0] != NULL) {
      free(user_elt_ln_to_gn_a[ipart]);
    }
  }
  free(user_elt_a_b_init_loc_idx);
  free(user_elt_a_b_init_loc    );
  free(user_elt_ln_to_gn_a);
  free(user_n_elt_a);
  free(user_n_elt_b);

  free(cellA_lineB_post);
  free(cellA_lineB_post_idx);
  free(cellA_lineB_post_n);

}



static void
_get_extracted_mesh_surf
(
 PDM_extract_part_t  *extrp,
 int                 *n_face,
 int                 *n_edge,
 int                 *n_vtx,
 int                **face_edge_idx,
 int                **face_edge,
 int                **edge_vtx,
 double             **vtx_coord,
 PDM_g_num_t        **face_ln_to_gn
 )
 {
  *n_face = PDM_extract_part_connectivity_get(extrp,
                                              0,
                                              PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                              face_edge,
                                              face_edge_idx,
                                              PDM_OWNERSHIP_KEEP);
  int *edge_vtx_idx = NULL;
  *n_edge = PDM_extract_part_connectivity_get(extrp,
                                              0,
                                              PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                              edge_vtx,
                                              &edge_vtx_idx,
                                              PDM_OWNERSHIP_KEEP);

  *n_vtx = PDM_extract_part_vtx_coord_get(extrp,
                                          0,
                                          vtx_coord,
                                          PDM_OWNERSHIP_KEEP);

  *face_ln_to_gn = extrp->target_gnum[0];
  // PDM_extract_part_parent_ln_to_gn_get(extrp,
  //                                      0,
  //                                      PDM_MESH_ENTITY_FACE,
  //                                      face_ln_to_gn,
  //                                      PDM_OWNERSHIP_KEEP);
 }


static inline void
_vector_ab
(
      double ab[3],
const double a[3],
const double b[3]
)
{
  ab[0] = b[0] - a[0];
  ab[1] = b[1] - a[1];
  ab[2] = b[2] - a[2];
}

static void
_polygon_geom_properties
(
 int     n_edge,
 int    *face_edge,
 int    *edge_vtx,
 double *vtx_coord,
 double *normal,
 double *barycenter
 )
{
  for (int i = 0; i < 3; i++) {
    normal    [i] = 0;
    barycenter[i] = 0;
  }

  for (int iedge = 0; iedge < n_edge; iedge++) {
    int edge_id   = PDM_ABS (face_edge[iedge]) - 1;

    int vtx_id0 = edge_vtx[2*edge_id  ] - 1;
    int vtx_id1 = edge_vtx[2*edge_id+1] - 1;

    for (int i = 0; i < 3; i++) {
      barycenter[i] += vtx_coord[3*vtx_id0+i] + vtx_coord[3*vtx_id1+i];
    }
  }

  double normalization = 1./(2. * n_edge);

  for (int i = 0; i < 3; i++) {
    barycenter[i] *= normalization;
  }

  for (int iedge = 0; iedge < n_edge; iedge++) {
    int edge_id   = PDM_ABS (face_edge[iedge]) - 1;
    int edge_sign = PDM_SIGN(face_edge[iedge]);

    double vec[2][3];
    for (int j = 0; j < 2; j++) {
      int vtx_id = edge_vtx[2*edge_id+j] - 1;
      _vector_ab(vec[j], barycenter, &vtx_coord[3*vtx_id]);
    }

    double cross[3];
    PDM_CROSS_PRODUCT(cross, vec[0], vec[1]);
    for (int i = 0; i < 3; i++) {
      normal[i] += edge_sign * cross[i];
    }
  }
}


PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
static inline void
_clip1
(
       double *uc,
       double *vc,
 const double  ud,
 const double  vd
 )
{
  if (*uc < 0) {
    if (*uc == ud) {
      *vc = 0;
    }
    else {
      *vc -= (*uc)*(vd - (*vc))/(ud - (*uc));
    }
    *uc = 0;
  }

  if (*vc < 0) {
    if (*vc == vd) {
      *uc = 0;
    }
    else {
      *uc -= (*vc)*(ud - (*uc))/(vd - (*vc));
    }
    *vc = 0;
  }
}

static inline void
_clip2
(
 double *u,
 double *v
 )
{
  double w = (*u) + (*v);
  if (w > 1) {
    double iw = 1./w;
    *u *= iw;
    *v *= iw;
  }
}
PDM_GCC_SUPPRESS_WARNING_POP



 static void
 _get_extracted_mesh_line
 (
  PDM_extract_part_t  *extrp,
  int                 *n_edge,
  int                 *n_vtx,
  int                **edge_vtx,
  double             **vtx_coord,
  PDM_g_num_t        **edge_ln_to_gn
  )
 {
  int *edge_vtx_idx = NULL;
  *n_edge = PDM_extract_part_connectivity_get(extrp,
                                              0,
                                              PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                              edge_vtx,
                                              &edge_vtx_idx,
                                              PDM_OWNERSHIP_KEEP);

  *n_vtx = PDM_extract_part_vtx_coord_get(extrp,
                                          0,
                                          vtx_coord,
                                          PDM_OWNERSHIP_KEEP);

  *edge_ln_to_gn = extrp->target_gnum[0];
  // PDM_extract_part_parent_ln_to_gn_get(extrp,
  //                                      0,
  //                                      PDM_MESH_ENTITY_EDGE,
  //                                      edge_ln_to_gn,
  //                                      PDM_OWNERSHIP_KEEP);
}


static
void
_mesh_intersection_surf_line
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  int i_rank;
  PDM_MPI_Comm_rank(mi->comm, &i_rank);

  int dbg_enabled = 0;

  if(dbg_enabled) {
    _export_vtk_2d("extrp_mesh_a", extrp_mesh_a);
    _export_vtk_1d("extrp_edge_a", extrp_mesh_a);
    _export_vtk_1d("extrp_mesh_b", extrp_mesh_b);
  }

  /* Get connectivities and coordinates */
  int          n_faceA         = 0;
  int          n_edgeA         = 0;
  int          n_vtxA          = 0;
  int         *faceA_edgeA_idx = NULL;
  int         *faceA_edgeA     = NULL;
  int         *edgeA_vtxA      = NULL;
  double      *vtxA_coord      = NULL;
  PDM_g_num_t *faceA_ln_to_gn  = NULL;
  _get_extracted_mesh_surf(extrp_mesh_a,
                           &n_faceA,
                           &n_edgeA,
                           &n_vtxA,
                           &faceA_edgeA_idx,
                           &faceA_edgeA,
                           &edgeA_vtxA,
                           &vtxA_coord,
                           &faceA_ln_to_gn);

  int          n_edgeB         = 0;
  int          n_vtxB          = 0;
  int         *edgeB_vtxB      = NULL;
  double      *vtxB_coord      = NULL;
  PDM_g_num_t *edgeB_ln_to_gn  = NULL;
  _get_extracted_mesh_line(extrp_mesh_b,
                           &n_edgeB,
                           &n_vtxB,
                           &edgeB_vtxB,
                           &vtxB_coord,
                           &edgeB_ln_to_gn);



  int *faceA_edgeB_idx = redistribute_box_a_to_box_b_idx;
  int *faceA_edgeB     = redistribute_box_a_to_box_b;


  /*
   We need :
     - subedgeB_vector/center/parent_edgeB/faceA
     - subedgeA_vector/center/parent_edgeA
  */

  int *edgeA_faceA = PDM_array_const_int(2*n_edgeA, -1);
  for (int faceA_id = 0; faceA_id < n_faceA; faceA_id++) {
    for (int iedgeA = faceA_edgeA_idx[faceA_id]; iedgeA < faceA_edgeA_idx[faceA_id+1]; iedgeA++) {
      int edgeA_id   = PDM_ABS (faceA_edgeA[iedgeA]) - 1;
      int edgeA_sign = PDM_SIGN(faceA_edgeA[iedgeA]);

      int idx = 2*edgeA_id;
      if (edgeA_sign < 0) {
        idx++;
      }

      edgeA_faceA[idx] = faceA_id;
    }
  }

  int *edgeA_inter_n = PDM_array_zeros_int(n_edgeA);
  int *edgeB_inter_n = PDM_array_zeros_int(n_edgeB);

  for (int faceA_id = 0; faceA_id < n_faceA; faceA_id++) {

    for (int iedgeB = faceA_edgeB_idx[faceA_id]; iedgeB < faceA_edgeB_idx[faceA_id+1]; iedgeB++) {

      int edgeB_id = faceA_edgeB[iedgeB];

      int vtxB_id0 = edgeB_vtxB[2*edgeB_id  ] - 1;
      int vtxB_id1 = edgeB_vtxB[2*edgeB_id+1] - 1;

      for (int iedgeA = faceA_edgeA_idx[faceA_id]; iedgeA < faceA_edgeA_idx[faceA_id+1]; iedgeA++) {

        int edgeA_id = PDM_ABS (faceA_edgeA[iedgeA]) - 1;

        int other_faceA = -1;
        for (int ifaceA = 2*edgeA_id; ifaceA < 2*(edgeA_id+1); ifaceA++) {
          if (edgeA_faceA[ifaceA] != faceA_id) {
            other_faceA = edgeA_faceA[ifaceA];
            break;
          }
        }

        if (other_faceA > faceA_id) {
          continue;
        }

        int vtxA_id0 = edgeA_vtxA[2*edgeA_id  ] - 1;
        int vtxA_id1 = edgeA_vtxA[2*edgeA_id+1] - 1;

        double tA, tB;
        int stat = PDM_line_intersection_mean_square(&vtxA_coord[3*vtxA_id0],
                                                     &vtxA_coord[3*vtxA_id1],
                                                     &vtxB_coord[3*vtxB_id0],
                                                     &vtxB_coord[3*vtxB_id1],
                                                     &tA,
                                                     &tB);

        if (stat == PDM_LINE_INTERSECT_YES) {
          edgeA_inter_n[edgeA_id]++;
          edgeB_inter_n[edgeB_id]++;
        }

      } // End of loop on edges A of current face A

    } // End of loop on candidate edges B for current face A

  } // End of loop on faces A

  int max_n = 0;
  for (int edgeA_id = 0; edgeA_id < n_edgeA; edgeA_id++) {
    max_n = PDM_MAX(max_n, edgeA_inter_n[edgeA_id]);
  }
  for (int edgeB_id = 0; edgeB_id < n_edgeB; edgeB_id++) {
    max_n = PDM_MAX(max_n, edgeB_inter_n[edgeB_id]);
  }


  int *edgeA_inter_idx = PDM_array_new_idx_from_sizes_int(edgeA_inter_n, n_edgeA);
  int *edgeB_inter_idx = PDM_array_new_idx_from_sizes_int(edgeB_inter_n, n_edgeB);

  PDM_array_reset_int(edgeA_inter_n, n_edgeA, 0);
  PDM_array_reset_int(edgeB_inter_n, n_edgeB, 0);


  typedef enum {
    ENTERING,
    EXITING
  } _crossing_t;


  double *edgeA_inter_t     = malloc(sizeof(double) * edgeA_inter_idx[n_edgeA]);
  int    *edgeA_inter_edgeB = malloc(sizeof(int   ) * edgeA_inter_idx[n_edgeA]);
  double *edgeB_inter_t     = malloc(sizeof(double) * edgeB_inter_idx[n_edgeB]);
  int    *edgeB_inter_edgeA = malloc(sizeof(int   ) * edgeB_inter_idx[n_edgeB]);

  _crossing_t *edgeA_inter_crossing = malloc(sizeof(int) * edgeA_inter_idx[n_edgeA]);
  _crossing_t *edgeB_inter_crossing = malloc(sizeof(int) * edgeB_inter_idx[n_edgeB]);

  for (int faceA_id = 0; faceA_id < n_faceA; faceA_id++) {

    for (int iedgeB = faceA_edgeB_idx[faceA_id]; iedgeB < faceA_edgeB_idx[faceA_id+1]; iedgeB++) {

      int edgeB_id = faceA_edgeB[iedgeB];

      int vtxB_id0 = edgeB_vtxB[2*edgeB_id  ] - 1;
      int vtxB_id1 = edgeB_vtxB[2*edgeB_id+1] - 1;

      for (int iedgeA = faceA_edgeA_idx[faceA_id]; iedgeA < faceA_edgeA_idx[faceA_id+1]; iedgeA++) {

        int edgeA_id = PDM_ABS (faceA_edgeA[iedgeA]) - 1;

        int other_faceA = -1;
        for (int ifaceA = 2*edgeA_id; ifaceA < 2*(edgeA_id+1); ifaceA++) {
          if (edgeA_faceA[ifaceA] != faceA_id) {
            other_faceA = edgeA_faceA[ifaceA];
            break;
          }
        }

        if (other_faceA > faceA_id) {
          continue;
        }

        int vtxA_id0 = edgeA_vtxA[2*edgeA_id  ] - 1;
        int vtxA_id1 = edgeA_vtxA[2*edgeA_id+1] - 1;

        double tA, tB;
        int stat = PDM_line_intersection_mean_square(&vtxA_coord[3*vtxA_id0],
                                                     &vtxA_coord[3*vtxA_id1],
                                                     &vtxB_coord[3*vtxB_id0],
                                                     &vtxB_coord[3*vtxB_id1],
                                                     &tA,
                                                     &tB);

        if (stat == PDM_LINE_INTERSECT_YES) {
          int idxA = edgeA_inter_idx[edgeA_id] + edgeA_inter_n[edgeA_id];
          int idxB = edgeB_inter_idx[edgeB_id] + edgeB_inter_n[edgeB_id];

          edgeA_inter_t[idxA] = tA;
          edgeB_inter_t[idxB] = tB;

          edgeA_inter_edgeB[idxA] = edgeB_id;
          edgeB_inter_edgeA[idxB] = edgeA_id;

          // /!\ if not in plane xy
          double detA = PDM_predicate_orient2d(&vtxB_coord[3*vtxB_id0],
                                               &vtxB_coord[3*vtxB_id1],
                                               &vtxA_coord[3*vtxA_id0]);
          if (detA < 0) {
            edgeA_inter_crossing[idxA] = ENTERING;
          }
          else {
            edgeA_inter_crossing[idxA] = EXITING;
          }

          // /!\ if not in plane xy
          double detB = PDM_predicate_orient2d(&vtxA_coord[3*vtxA_id0],
                                               &vtxA_coord[3*vtxA_id1],
                                               &vtxB_coord[3*vtxB_id0]);
          if (detB < 0) {
            edgeB_inter_crossing[idxB] = ENTERING;
          }
          else {
            edgeB_inter_crossing[idxB] = EXITING;
          }

          edgeA_inter_n[edgeA_id]++;
          edgeB_inter_n[edgeB_id]++;
        }

      } // End of loop on edges A of current face A

    } // End of loop on candidate edges B for current face A

  } // End of loop on faces A



  /* Sort intersection points along each edge */
  int *order = malloc(sizeof(int) * max_n);

  int s_subedgeA = n_edgeA + edgeA_inter_idx[n_edgeA];
  int    *subedgeA_parent = malloc(sizeof(int   ) * s_subedgeA);
  double *subedgeA_center = malloc(sizeof(double) * s_subedgeA * 3);
  double *subedgeA_vector = malloc(sizeof(double) * s_subedgeA * 3);


  double *dbg_subedgeA_coord = NULL;
  if (dbg_enabled) {
    dbg_subedgeA_coord = malloc(sizeof(double) * s_subedgeA * 6);
  }


  int n_subedgeA = 0;

  for (int edgeA_id = 0; edgeA_id < n_edgeA; edgeA_id++) {
    int vtxA_id0 = edgeA_vtxA[2*edgeA_id  ] - 1;
    int vtxA_id1 = edgeA_vtxA[2*edgeA_id+1] - 1;

    if (dbg_enabled) {
      log_trace("edgeA %d (%d %d):\n", edgeA_id, vtxA_id0, vtxA_id1);
    }

    int n = edgeA_inter_idx[edgeA_id+1] - edgeA_inter_idx[edgeA_id];
    int         *_edgeA_inter_edgeB    = edgeA_inter_edgeB    + edgeA_inter_idx[edgeA_id];
    double      *_edgeA_inter_t        = edgeA_inter_t        + edgeA_inter_idx[edgeA_id];
    _crossing_t *_edgeA_inter_crossing = edgeA_inter_crossing + edgeA_inter_idx[edgeA_id];

    if (n > 0) {

      for (int i = 0; i < n; i++) {
        order[i] = i;
      }

      PDM_sort_double(_edgeA_inter_t, order, n);

      int start = 0;
      if (_edgeA_inter_crossing[order[0]] == EXITING) {
        start++;
      }


      double p0[3];
      double p1[3];
      for (int i = start; i <= n; i+=2) {

        if (i == 0) {
          memcpy(p0, &vtxA_coord[3*vtxA_id0], sizeof(double)*3);
        }
        else {
          for (int j = 0; j < 3; j++) {
            p0[j] = (1-_edgeA_inter_t[i-1])*vtxA_coord[3*vtxA_id0+j] + _edgeA_inter_t[i-1]*vtxA_coord[3*vtxA_id1+j];
          }
        }

        if (i == n) {
          memcpy(p1, &vtxA_coord[3*vtxA_id1], sizeof(double)*3);
        }
        else {
          if (i < n-1) {
            assert(_edgeA_inter_crossing[order[i+1]] == ENTERING);
          }
          for (int j = 0; j < 3; j++) {
            p1[j] = (1-_edgeA_inter_t[i])*vtxA_coord[3*vtxA_id0+j] + _edgeA_inter_t[i]*vtxA_coord[3*vtxA_id1+j];
          }
        }

        subedgeA_parent[n_subedgeA] = edgeA_id;
        for (int j = 0; j < 3; j++) {
          subedgeA_center[3*n_subedgeA+j] = 0.5*(p0[j] + p1[j]);
        }
        // /!\ if not in plane xy
        subedgeA_vector[3*n_subedgeA  ] = p1[1] - p0[1];
        subedgeA_vector[3*n_subedgeA+1] = p0[0] - p1[0];
        subedgeA_vector[3*n_subedgeA+2] = 0;

        if (dbg_enabled) {
          memcpy(&dbg_subedgeA_coord[6*n_subedgeA  ], p0, sizeof(double)*3);
          memcpy(&dbg_subedgeA_coord[6*n_subedgeA+3], p1, sizeof(double)*3);
        }

        n_subedgeA++;
      }

    }

    if (dbg_enabled) {

      for (int i = 0; i < n; i++) {
        int edgeB_id = _edgeA_inter_edgeB[order[i]];
        log_trace("  edgeB %d (%d %d), t = %f, crossing %d\n",
                  edgeB_id,
                  edgeB_vtxB[2*edgeB_id], edgeB_vtxB[2*edgeB_id+1],
                  _edgeA_inter_t[i],
                  (int) _edgeA_inter_crossing[order[i]]);
      }
    }
  }

  if (dbg_enabled) {
    char filename[999];
    sprintf(filename, "dbg_subedgeA_%d.vtk", i_rank);

    PDM_vtk_write_lines(filename,
                        n_subedgeA,
                        dbg_subedgeA_coord,
                        NULL,
                        subedgeA_parent);

    sprintf(filename, "dbg_subedgeA_%d_vector.vtk", i_rank);
    const char   *vector_field_name[1] = {"vector"};
    const double *vector_field     [1] = {subedgeA_vector};
    PDM_vtk_write_point_cloud_with_field(filename,
                                         n_subedgeA,
                                         subedgeA_center,
                                         NULL,
                                         subedgeA_parent,
                                         0,
                                         NULL,
                                         NULL,
                                         1,
                                         vector_field_name,
                                         vector_field,
                                         0,
                                         NULL,
                                         NULL);

    free(dbg_subedgeA_coord);
  }

  if (dbg_enabled) {
    log_trace("\n\n");
  }



  int s_subedgeB = n_edgeB + edgeB_inter_idx[n_edgeB];
  int    *subedgeB_parent = malloc(sizeof(int   ) * s_subedgeB);
  int    *subedgeB_faceA  = malloc(sizeof(int   ) * s_subedgeB);
  double *subedgeB_center = malloc(sizeof(double) * s_subedgeB * 3);
  double *subedgeB_vector = malloc(sizeof(double) * s_subedgeB * 3);


  double *dbg_subedgeB_coord = NULL;
  if (dbg_enabled) {
    dbg_subedgeB_coord = malloc(sizeof(double) * s_subedgeB * 6);
  }


  int n_subedgeB = 0;
  int n_undefB   = 0;

  for (int edgeB_id = 0; edgeB_id < n_edgeB; edgeB_id++) {
    int vtxB_id0 = edgeB_vtxB[2*edgeB_id  ] - 1;
    int vtxB_id1 = edgeB_vtxB[2*edgeB_id+1] - 1;

    if (dbg_enabled) {
      log_trace("edgeB %d (%d %d):\n", edgeB_id, vtxB_id0, vtxB_id1);
    }

    int n = edgeB_inter_idx[edgeB_id+1] - edgeB_inter_idx[edgeB_id];
    int         *_edgeB_inter_edgeA    = edgeB_inter_edgeA    + edgeB_inter_idx[edgeB_id];
    double      *_edgeB_inter_t        = edgeB_inter_t        + edgeB_inter_idx[edgeB_id];
    _crossing_t *_edgeB_inter_crossing = edgeB_inter_crossing + edgeB_inter_idx[edgeB_id];

    if (n == 0) {
      /* Current edge B is either outside A or inside a single face A */
      // we will deal with it later...
      n_undefB++;
    }
    else {
      for (int i = 0; i < n; i++) {
        order[i] = i;
      }

      PDM_sort_double(_edgeB_inter_t, order, n);

      double p0[3];
      double p1[3];
      int faceA_id;
      for (int i = 0; i <= n; i++) {

        if (i == 0) {
          memcpy(p0, &vtxB_coord[3*vtxB_id0], sizeof(double)*3);

          int edgeA_id = _edgeB_inter_edgeA[order[0]];
          if (_edgeB_inter_crossing[order[0]] == ENTERING) {
            faceA_id = edgeA_faceA[2*edgeA_id+1];
          }
          else {
            faceA_id = edgeA_faceA[2*edgeA_id  ];
          }
        }
        else {
          for (int j = 0; j < 3; j++) {
            p0[j] = (1-_edgeB_inter_t[i-1])*vtxB_coord[3*vtxB_id0+j] + _edgeB_inter_t[i-1]*vtxB_coord[3*vtxB_id1+j];
          }

          int edgeA_id = _edgeB_inter_edgeA[order[i-1]];
          if (_edgeB_inter_crossing[order[i-1]] == EXITING) {
            faceA_id = edgeA_faceA[2*edgeA_id+1];
          }
          else {
            faceA_id = edgeA_faceA[2*edgeA_id  ];
          }
        }

        if (i == n) {
          memcpy(p1, &vtxB_coord[3*vtxB_id1], sizeof(double)*3);
        }
        else {
          for (int j = 0; j < 3; j++) {
            p1[j] = (1-_edgeB_inter_t[i])*vtxB_coord[3*vtxB_id0+j] + _edgeB_inter_t[i]*vtxB_coord[3*vtxB_id1+j];
          }
        }

        if (faceA_id >= 0) {

          subedgeB_parent[n_subedgeB] = edgeB_id;
          subedgeB_faceA [n_subedgeB] = faceA_id;
          for (int j = 0; j < 3; j++) {
            subedgeB_center[3*n_subedgeB+j] = 0.5*(p0[j] + p1[j]);
          }
          // /!\ if not in plane xy
          subedgeB_vector[3*n_subedgeB  ] = p1[1] - p0[1];
          subedgeB_vector[3*n_subedgeB+1] = p0[0] - p1[0];
          subedgeB_vector[3*n_subedgeB+2] = 0;

          if (dbg_enabled) {
            log_trace("+ subedge %d : faceA_id = %d\n", n_subedgeB, faceA_id);
          }

          if (dbg_enabled) {
            memcpy(&dbg_subedgeB_coord[6*n_subedgeB  ], p0, sizeof(double)*3);
            memcpy(&dbg_subedgeB_coord[6*n_subedgeB+3], p1, sizeof(double)*3);
          }

          n_subedgeB++;
        }
      }
    }

    if (dbg_enabled) {
      for (int i = 0; i < n; i++) {
        int edgeA_id = _edgeB_inter_edgeA[order[i]];
        log_trace("  edgeA %d (%d %d), t = %f, crossing %d\n",
                  edgeA_id,
                  edgeA_vtxA[2*edgeA_id], edgeA_vtxA[2*edgeA_id+1],
                  _edgeB_inter_t[i],
                  (int) _edgeB_inter_crossing[order[i]]);
      }
    }
  }
  free(order);

  free(edgeA_inter_n);
  free(edgeA_inter_idx);
  free(edgeA_inter_t);
  free(edgeA_inter_edgeB);
  free(edgeA_inter_crossing);


  /* Deal with 'undef' edges B */
  if (n_undefB > 0) {

    /* Localize the undef edgeB midpoints */
    int *edgeB_faceA = PDM_array_zeros_int(n_edgeB);
    for (int edgeB_id = 0; edgeB_id < n_edgeB; edgeB_id++) {

      int n = edgeB_inter_idx[edgeB_id+1] - edgeB_inter_idx[edgeB_id];

      if (n > 0) {
        edgeB_faceA[edgeB_id] = -1;
      }
    }

    int *faceA_vtxA = NULL;
    PDM_compute_face_vtx_from_face_and_edge(n_faceA,
                                            faceA_edgeA_idx,
                                            faceA_edgeA,
                                            edgeA_vtxA,
                                            &faceA_vtxA);

    int max_face_vtx_n = 0;
    for (int faceA_id = 0; faceA_id < n_faceA; faceA_id++) {
      int n = faceA_edgeA_idx[faceA_id+1] - faceA_edgeA_idx[faceA_id];
      max_face_vtx_n = PDM_MAX(max_face_vtx_n, n);
    }


    double *faceA_coord = malloc(sizeof(double) * max_face_vtx_n * 3);

    for (int faceA_id = 0; faceA_id < n_faceA; faceA_id++) {

      int faceA_vtxA_n = faceA_edgeA_idx[faceA_id+1] - faceA_edgeA_idx[faceA_id];
      for (int ivtxA = 0; ivtxA < faceA_vtxA_n; ivtxA++) {
        int vtxA_id = faceA_vtxA[faceA_edgeA_idx[faceA_id] + ivtxA] - 1;
        memcpy(&faceA_coord[3*ivtxA], &vtxA_coord[3*vtxA_id], sizeof(double)*3);
      }

      double faceA_normal[3];
      PDM_plane_normal(faceA_vtxA_n, faceA_coord, faceA_normal);

      for (int iedgeB = faceA_edgeB_idx[faceA_id]; iedgeB < faceA_edgeB_idx[faceA_id+1]; iedgeB++) {

        int edgeB_id = faceA_edgeB[iedgeB];

        if (edgeB_faceA[edgeB_id] != 0) {
          continue;
        }

        int vtxB_id0 = edgeB_vtxB[2*edgeB_id  ] - 1;
        int vtxB_id1 = edgeB_vtxB[2*edgeB_id+1] - 1;

        double edgeB_center[3];
        for (int i = 0; i < 3; i++) {
          edgeB_center[i] = 0.5*(vtxB_coord[3*vtxB_id0+i] + vtxB_coord[3*vtxB_id1+i]);
        }

        PDM_polygon_status_t stat = PDM_polygon_point_in_new(edgeB_center,
                                                             faceA_vtxA_n,
                                                             faceA_coord,
                                                             NULL,
                                                             faceA_normal);

        if (stat == PDM_POLYGON_INSIDE) {
          edgeB_faceA[edgeB_id] = faceA_id+1;

          subedgeB_parent[n_subedgeB] = edgeB_id;
          subedgeB_faceA [n_subedgeB] = faceA_id;
          memcpy(&subedgeB_center[3*n_subedgeB], edgeB_center, sizeof(double)*3);

          // /!\ if not in plane xy
          subedgeB_vector[3*n_subedgeB  ] = vtxB_coord[3*vtxB_id1+1] - vtxB_coord[3*vtxB_id0+1];
          subedgeB_vector[3*n_subedgeB+1] = vtxB_coord[3*vtxB_id0+0] - vtxB_coord[3*vtxB_id1+0];
          subedgeB_vector[3*n_subedgeB+2] = 0;

          if (dbg_enabled) {
            log_trace("+ subedge %d : faceA_id = %d\n", n_subedgeB, faceA_id);
          }

          if (dbg_enabled) {
            memcpy(&dbg_subedgeB_coord[6*n_subedgeB  ], &vtxB_coord[3*vtxB_id0], sizeof(double)*3);
            memcpy(&dbg_subedgeB_coord[6*n_subedgeB+3], &vtxB_coord[3*vtxB_id1], sizeof(double)*3);
          }

          n_subedgeB++;

        }

      }

    }

    free(edgeB_faceA);
    free(faceA_vtxA);
    free(faceA_coord);
  }

  free(edgeB_inter_n);
  free(edgeB_inter_idx);
  free(edgeB_inter_t);
  free(edgeB_inter_edgeA);
  free(edgeB_inter_crossing);



  if (dbg_enabled) {
    char filename[999];
    sprintf(filename, "dbg_subedgeB_%d.vtk", i_rank);

    PDM_vtk_write_lines(filename,
                        n_subedgeB,
                        dbg_subedgeB_coord,
                        NULL,
                        subedgeB_faceA);//subedgeB_parent);

    sprintf(filename, "dbg_subedgeB_%d_vector.vtk", i_rank);
    const char   *vector_field_name[1] = {"vector"};
    const double *vector_field     [1] = {subedgeB_vector};
    PDM_vtk_write_point_cloud_with_field(filename,
                                         n_subedgeB,
                                         subedgeB_center,
                                         NULL,
                                         subedgeB_faceA,
                                         0,
                                         NULL,
                                         NULL,
                                         1,
                                         vector_field_name,
                                         vector_field,
                                         0,
                                         NULL,
                                         NULL);

    free(dbg_subedgeB_coord);
  }






  if (dbg_enabled) {
    /* Check surface bilan */
    double *faceA_bilan = malloc(sizeof(double) * n_faceA * 3);
    for (int i = 0; i < 3*n_faceA; i++) {
      faceA_bilan[i] = 0;
    }

    for (int i = 0; i < n_subedgeA; i++) {
      int edgeA_id = subedgeA_parent[i];

      int sign = 1;
      for (int ifaceA = 2*edgeA_id; ifaceA < 2*(edgeA_id+1); ifaceA++) {
        int faceA_id = edgeA_faceA[ifaceA];

        if (faceA_id >= 0) {
          for (int j = 0; j < 3; j++) {
            faceA_bilan[3*faceA_id+j] += sign*subedgeA_vector[3*i+j];
          }
        }

        sign = -sign;
      }
    }


    for (int i = 0; i < n_subedgeB; i++) {
      int faceA_id = subedgeB_faceA[i];
      for (int j = 0; j < 3; j++) {
        faceA_bilan[3*faceA_id+j] -= subedgeB_vector[3*i+j];
      }
    }


    for (int i = 0; i < n_faceA; i++) {
      double mag = PDM_MODULE(&faceA_bilan[3*i]);
      log_trace("faceA %d : bilan = %e\n", i, mag);
    }

    free(faceA_bilan);
  }


  free(edgeA_faceA);




  /* Results, do not free */
  free(subedgeA_parent);
  free(subedgeA_center);
  free(subedgeA_vector);

  free(subedgeB_parent);
  free(subedgeB_faceA );
  free(subedgeB_center);
  free(subedgeB_vector);


}





static
void
_mesh_intersection_surf_surf
(
 PDM_mesh_intersection_t *mi,
 int                     *a_to_b_idx,
 int                     *a_to_b
)
{
  int i_rank;
  PDM_MPI_Comm_rank(mi->comm, &i_rank);

  int dbg_enabled = 0;

  PDM_part_mesh_nodal_elmts_t *pmne[2] = {NULL, NULL};



  int          n_vtx        [2];
  double      *vtx_coord    [2] = {NULL, NULL};
  int          n_face       [2];
  int         *face_vtx_idx [2] = {NULL, NULL};
  int         *face_vtx     [2] = {NULL, NULL};
  int         *face_edge    [2] = {NULL, NULL};
  int         *edge_vtx     [2] = {NULL, NULL};
  PDM_g_num_t *face_ln_to_gn[2] = {NULL, NULL};

  /* Get extract surface meshes */
  for (int i = 0; i < 2; i++) {

    n_vtx[i] = PDM_extract_part_n_entity_get(mi->extrp_mesh[i],
                                             0,
                                             PDM_MESH_ENTITY_VTX);

    PDM_extract_part_vtx_coord_get(mi->extrp_mesh[i],
                                   0,
                                   &vtx_coord[i],
                                   PDM_OWNERSHIP_KEEP);

    n_face[i] = mi->extrp_mesh[i]->n_target[0];
    face_ln_to_gn[i] = mi->extrp_mesh[i]->target_gnum[0];

    if (dbg_enabled) {
      log_trace("mesh %d is nodal? %d\n", i, mi->mesh[i] == NULL);
    }

    if (mi->mesh[i] == NULL) {
      PDM_extract_part_part_mesh_nodal_get(mi->extrp_mesh[i],
                                           &pmne[i],
                                           PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_nodal_elmts_t *extract_pmne = pmne[i];


      int n_section = PDM_part_mesh_nodal_elmts_n_section_get(extract_pmne);
      int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(extract_pmne);

      // n_face[i] = 0;
      // for (int isection = 0; isection < n_section; isection++) {
      //   n_face[i] += PDM_part_mesh_nodal_elmts_section_n_elt_get(extract_pmne,
      //                                                          sections_id[isection],
      //                                                          0);
      // }

      face_vtx_idx[i] = malloc(sizeof(int) * (n_face[i] + 1));
      face_vtx_idx[i][0] = 0;

      for (int isection = 0; isection < n_section; isection++) {
        int id_section = sections_id[isection];

        int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(extract_pmne,
                                                              id_section,
                                                              0);

        PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extract_pmne,
                                                                              id_section);

        assert(PDM_Mesh_nodal_elt_dim_get(t_elt) == 2);

        int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(extract_pmne,
                                                                   id_section,
                                                                   0,
                                                                   PDM_OWNERSHIP_KEEP);
        if (t_elt == PDM_MESH_NODAL_POLY_2D) {
          /* Polygonal section */
          int *connec_idx;
          int *connec;
          PDM_part_mesh_nodal_elmts_section_poly2d_get(extract_pmne,
                                                     id_section,
                                                     0,
                                                     &connec_idx,
                                                     &connec,
                                                     PDM_OWNERSHIP_KEEP);

          for (int ielt = 0; ielt < n_elt; ielt++) {
            int iface = ielt;
            if (parent_num != NULL) {
              iface = parent_num[ielt];
            }

            int face_vtx_n = connec_idx[ielt+1] - connec_idx[ielt];

            face_vtx_idx[i][iface+1] = face_vtx_n;
          }
        }
        else {
          /* Standard section */
          int         *connec              = NULL;
          PDM_g_num_t *numabs              = NULL;
          int         *_parent_num         = NULL;
          PDM_g_num_t *parent_entity_g_num = NULL;
          int          order               = 0;
          const char  *ho_ordering         = NULL;
          PDM_part_mesh_nodal_elmts_section_std_ho_get(extract_pmne,
                                                     id_section,
                                                     0,
                                                     &connec,
                                                     &numabs,
                                                     &_parent_num,
                                                     &parent_entity_g_num,
                                                     &order,
                                                     &ho_ordering,
                                                     PDM_OWNERSHIP_KEEP);
          assert(order == 1);

          int face_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                                        order);


          for (int ielt = 0; ielt < n_elt; ielt++) {
            int iface = ielt;
            if (parent_num != NULL) {
              iface = parent_num[ielt];
            }

            face_vtx_idx[i][iface+1] = face_vtx_n;
          }
        }

      } // End of loop on sections

      for (int iface = 0; iface < n_face[i]; iface++) {
        face_vtx_idx[i][iface+1] += face_vtx_idx[i][iface];
      }

      /* Fill face_vtx */
      face_vtx[i] = malloc(sizeof(int) * face_vtx_idx[i][n_face[i]]);
      for (int isection = 0; isection < n_section; isection++) {
        int id_section = sections_id[isection];

        int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(extract_pmne,
                                                              id_section,
                                                              0);

        PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(extract_pmne,
                                                                              id_section);

        assert(PDM_Mesh_nodal_elt_dim_get(t_elt) == 2);

        int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(extract_pmne,
                                                                   id_section,
                                                                   0,
                                                                   PDM_OWNERSHIP_KEEP);
        if (t_elt == PDM_MESH_NODAL_POLY_2D) {
          /* Polygonal section */
          int *connec_idx;
          int *connec;
          PDM_part_mesh_nodal_elmts_section_poly2d_get(extract_pmne,
                                                     id_section,
                                                     0,
                                                     &connec_idx,
                                                     &connec,
                                                     PDM_OWNERSHIP_KEEP);

          for (int ielt = 0; ielt < n_elt; ielt++) {
            int iface = ielt;
            if (parent_num != NULL) {
              iface = parent_num[ielt];
            }

            int face_vtx_n = connec_idx[ielt+1] - connec_idx[ielt];
            for (int j = 0; j < face_vtx_n; j++) {
              face_vtx[i][face_vtx_idx[i][iface]+j] = connec[connec_idx[ielt]+j];
            }
          }
        }
        else {
          /* Standard section */
          int         *connec              = NULL;
          PDM_g_num_t *numabs              = NULL;
          int         *_parent_num         = NULL;
          PDM_g_num_t *parent_entity_g_num = NULL;
          int          order               = 0;
          const char  *ho_ordering         = NULL;
          PDM_part_mesh_nodal_elmts_section_std_ho_get(extract_pmne,
                                                     id_section,
                                                     0,
                                                     &connec,
                                                     &numabs,
                                                     &_parent_num,
                                                     &parent_entity_g_num,
                                                     &order,
                                                     &ho_ordering,
                                                     PDM_OWNERSHIP_KEEP);
          assert(order == 1);

          int face_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(t_elt,
                                                        order);


          for (int ielt = 0; ielt < n_elt; ielt++) {
            int iface = ielt;
            if (parent_num != NULL) {
              iface = parent_num[ielt];
            }

            for (int j = 0; j < face_vtx_n; j++) {
              face_vtx[i][face_vtx_idx[i][iface]+j] = connec[face_vtx_n*ielt+j];
            }
          }
        }

      } // End of loop on sections
    }
    else {
      PDM_extract_part_connectivity_get(mi->extrp_mesh[i],
                                        0,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &face_vtx[i],
                                        &face_vtx_idx[i],
                                        PDM_OWNERSHIP_KEEP);
      int *_face_edge_idx = NULL;
      PDM_extract_part_connectivity_get(mi->extrp_mesh[i],
                                        0,
                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                        &face_edge[i],
                                        &_face_edge_idx,
                                        PDM_OWNERSHIP_KEEP);
      if (face_vtx_idx[i] == NULL) {
        assert(_face_edge_idx != NULL);
        face_vtx_idx[i] = _face_edge_idx;
      }

      int *_edge_vtx_idx = NULL;
      PDM_extract_part_connectivity_get(mi->extrp_mesh[i],
                                        0,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx[i],
                                        &_edge_vtx_idx,
                                        PDM_OWNERSHIP_KEEP);


    }

    if (dbg_enabled) {
      char filename[999];
      sprintf(filename, "mesh_intersection_surf_surf_mesh%d_rank%d.vtk", i, i_rank);

      int *_face_vtx = face_vtx[i];
      if (face_vtx[i] == NULL) {
        PDM_compute_face_vtx_from_face_and_edge(n_face[i],
                                                face_vtx_idx[i],
                                                face_edge[i],
                                                edge_vtx[i],
                                                &_face_vtx);
      }

      PDM_vtk_write_polydata(filename,
                             n_vtx[i],
                             vtx_coord[i],
                             NULL,
                             n_face[i],
                             face_vtx_idx[i],
                             _face_vtx,
                             face_ln_to_gn[i],
                             NULL);

      if (face_vtx[i] == NULL) {
        free(_face_vtx);
      }
    }
  }


  double *a_to_b_volume = malloc(sizeof(double) * a_to_b_idx[n_face[0]]);


  /* Compute face normals of B */
  double *faceB_normals = malloc(sizeof(double) * n_face[1] * 3);
  for (int faceB_id = 0; faceB_id < n_face[1]; faceB_id++) {
    double faceB_center[3];
    if (face_vtx[1] == NULL) {
      _polygon_geom_properties(face_vtx_idx[1][faceB_id+1] - face_vtx_idx[1][faceB_id],
                               face_edge[1] + face_vtx_idx[1][faceB_id],
                               edge_vtx[1],
                               vtx_coord[1],
                               faceB_normals + 3*faceB_id,
                               faceB_center);
    }
    else {
      int is_degenerate = 0;
      PDM_geom_elem_polygon_properties(1,
                                       face_vtx_idx[1] + faceB_id,
                                       face_vtx[1],
                                       vtx_coord[1],
                                       faceB_normals + 3*faceB_id,
                                       faceB_center,
                                       NULL,
                                       &is_degenerate);
      if (is_degenerate) {
        // PDM_error()
      }
    }
  }

  /* Main loop */
  for (int faceA_id = 0; faceA_id < n_face[0]; faceA_id++) {

    /* Compute faceA center and normal vector */
    double faceA_normal[3];
    double faceA_center[3];
    if (face_vtx[0] == NULL) {
      _polygon_geom_properties(face_vtx_idx[0][faceA_id+1] - face_vtx_idx[0][faceA_id],
                               face_edge[0] + face_vtx_idx[0][faceA_id],
                               edge_vtx[0],
                               vtx_coord[0],
                               faceA_normal,
                               faceA_center);
    }
    else {
      int is_degenerate = 0;
      PDM_geom_elem_polygon_properties(1,
                                       face_vtx_idx[0] + faceA_id,
                                       face_vtx[0],
                                       vtx_coord[0],
                                       faceA_normal,
                                       faceA_center,
                                       NULL,
                                       &is_degenerate);

      if (is_degenerate) {
        continue;
      }
    }


    /* Unit normal */
    double mag = PDM_DOT_PRODUCT(faceA_normal, faceA_normal);
    if (mag <= 0) {
      // degenerate polygon
      continue;
    }

    double imag = 1./sqrt(mag);
    for (int i = 0; i < 3; i++) {
      faceA_normal[i] *= imag;
    }

    int faceA_vtx_n = face_vtx_idx[0][faceA_id+1] - face_vtx_idx[0][faceA_id];


    for (int ifaceB = a_to_b_idx[faceA_id]; ifaceB < a_to_b_idx[faceA_id+1]; ifaceB++) {

      int faceB_id = a_to_b[ifaceB];

      int dbg_pair = 0;// dbg_enabled && (faceA_id == 5 && faceB_id == 0);

      double area = 0.;

      int signAB = (int) PDM_SIGN(PDM_DOT_PRODUCT(faceA_normal, faceB_normals + 3*faceB_id));
      if (dbg_pair) {
        log_trace("faceA %d ("PDM_FMT_G_NUM") faceB %d ("PDM_FMT_G_NUM"), signAB = %d\n",
                  faceA_id, face_ln_to_gn[0][faceA_id],
                  faceB_id, face_ln_to_gn[1][faceB_id],
                  signAB);
        log_trace("faceA_center = %f %f %f\n", faceA_center[0], faceA_center[1], faceA_center[2]);
        log_trace("faceA_normal = %f %f %f\n", faceA_normal[0], faceA_normal[1], faceA_normal[2]);
      }

      int faceB_vtx_n = face_vtx_idx[1][faceB_id+1] - face_vtx_idx[1][faceB_id];


      for (int idxA = 0; idxA < faceA_vtx_n; idxA++) {

        int vtxA_id0 = -1;
        int vtxA_id1 = -1;

        if (face_vtx[0] == NULL) {
          int iedgeA = face_vtx_idx[0][faceA_id] + idxA;
          int edgeA_id   = PDM_ABS (face_edge[0][iedgeA]) - 1;
          int edgeA_sign = PDM_SIGN(face_edge[0][iedgeA]);

          vtxA_id0 = edge_vtx[0][2*edgeA_id  ] - 1;
          vtxA_id1 = edge_vtx[0][2*edgeA_id+1] - 1;
          if (edgeA_sign < 0) {
            int tmp = vtxA_id0;
            vtxA_id0 = vtxA_id1;
            vtxA_id1 = tmp;
          }
        }
        else {
          vtxA_id0 = face_vtx[0][face_vtx_idx[0][faceA_id] + idxA]                 - 1;
          vtxA_id1 = face_vtx[0][face_vtx_idx[0][faceA_id] + (idxA+1)%faceA_vtx_n] - 1;
        }

        double *a = vtx_coord[0] + 3*vtxA_id0;
        double *b = vtx_coord[0] + 3*vtxA_id1;


        double ka[3], kb[3];
        _vector_ab(ka, faceA_center, a);
        _vector_ab(kb, faceA_center, b);

        double kaka = PDM_DOT_PRODUCT(ka, ka);
        double kakb = PDM_DOT_PRODUCT(ka, kb);
        double kbkb = PDM_DOT_PRODUCT(kb, kb);

        double det = kaka*kbkb - kakb*kakb;

        if (det <= 0.) {
          // points k, a and b are collinear, skip edge ab
          continue;
        }

        double idet = 1./det;

        double normal_kab[3];
        PDM_CROSS_PRODUCT(normal_kab, ka, kb);

        double area_kab = 0.5*PDM_DOT_PRODUCT(normal_kab, faceA_normal);

        for (int idxB = 0; idxB < faceB_vtx_n; idxB++) {

          int vtxB_id0 = -1;
          int vtxB_id1 = -1;

          if (face_vtx[1] == NULL) {
            int iedgeB = face_vtx_idx[1][faceB_id] + idxB;
            int edgeB_id   = PDM_ABS (face_edge[1][iedgeB]) - 1;
            int edgeB_sign = PDM_SIGN(face_edge[1][iedgeB]);

            vtxB_id0 = edge_vtx[1][2*edgeB_id  ] - 1;
            vtxB_id1 = edge_vtx[1][2*edgeB_id+1] - 1;
            if (edgeB_sign < 0) {
              int tmp = vtxB_id0;
              vtxB_id0 = vtxB_id1;
              vtxB_id1 = tmp;
            }
          }
          else {
            vtxB_id0 = face_vtx[1][face_vtx_idx[1][faceB_id] + idxB]                 - 1;
            vtxB_id1 = face_vtx[1][face_vtx_idx[1][faceB_id] + (idxB+1)%faceB_vtx_n] - 1;
          }

          double *c = vtx_coord[1] + 3*vtxB_id0;
          double *d = vtx_coord[1] + 3*vtxB_id1;

          /* Compute the barycentric coordinates of c and d in kab */
          double kc[3], kd[3];
          _vector_ab(kc, faceA_center, c);
          _vector_ab(kd, faceA_center, d);

          double kcka = PDM_DOT_PRODUCT(kc, ka);
          double kckb = PDM_DOT_PRODUCT(kc, kb);
          double kdka = PDM_DOT_PRODUCT(kd, ka);
          double kdkb = PDM_DOT_PRODUCT(kd, kb);

          double uc = kcka*kbkb - kckb*kakb;
          double ud = kdka*kbkb - kdkb*kakb;

          if (dbg_pair) {
            log_trace("      uc = %f, ud = %f\n", uc, ud);
          }

          if (uc <= 0 && ud <= 0) {
            continue;
          }

          double vc = kaka*kckb - kakb*kcka;
          double vd = kaka*kdkb - kakb*kdka;

          if (dbg_pair) {
            log_trace("      vc = %f, vd = %f\n", vc, vd);
          }

          if (vc <= 0 && vd <= 0) {
            continue;
          }

          PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
          if (uc*vd - vc*ud == 0) {
            // points k, c and d are collinear, skip edge cd
            continue;
          }
          PDM_GCC_SUPPRESS_WARNING_POP

          uc *= idet;
          vc *= idet;
          ud *= idet;
          vd *= idet;

          /* Compute intersection between ab and cd */
          int intersect = 0;
          double m00 = -1;
          double m01 = uc - ud;
          double m10 = 1;
          double m11 = vc - vd;
          double s;

          double det2 = m00*m11 - m01*m10;
          PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
          if (det2 != 0) {
            double idet2 = 1./det2;
            double r0 = uc - 1;
            double r1 = vc;
            s        = (r0*m11 - r1*m01) * idet2;
            double t = (m00*r1 - m10*r0) * idet2;
            if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
              intersect = 1;
            }
          }
          PDM_GCC_SUPPRESS_WARNING_POP

          /* Clip kcd by 'quarter space' {u,v >= 0} */
          _clip1(&uc, &vc, ud, vd);
          _clip1(&ud, &vd, uc, vc);

          /* Clip kcd by triangle kab {u+v <= 1} */
          _clip2(&uc, &vc);
          _clip2(&ud, &vd);

          if (dbg_pair) {
            log_trace("      final uc = %f, vc = %f\n", uc, vc);
            log_trace("      final ud = %f, vd = %f\n", ud, vd);
          }

          /* Add contribution */
          double f = 0;
          if (intersect) {
            f = (1-s)*(vd - vc) + s*(uc - ud);
            if (dbg_pair) {
              log_trace("    c2 = %f %f %f\n",
                        faceA_center[0] + uc*ka[0] + vc*kb[0],
                        faceA_center[1] + uc*ka[1] + vc*kb[1],
                        faceA_center[2] + uc*ka[2] + vc*kb[2]);
              log_trace("    x  = %f %f %f\n",
                        faceA_center[0] + (1-s)*ka[0] + s*kb[0],
                        faceA_center[1] + (1-s)*ka[1] + s*kb[1],
                        faceA_center[2] + (1-s)*ka[2] + s*kb[2]);
              log_trace("    d2 = %f %f %f\n",
                        faceA_center[0] + ud*ka[0] + vd*kb[0],
                        faceA_center[1] + ud*ka[1] + vd*kb[1],
                        faceA_center[2] + ud*ka[2] + vd*kb[2]);
            }
          }
          else {
            f = uc*vd - vc*ud;
            if (dbg_pair) {
              log_trace("    c2 = %f %f %f\n",
                        faceA_center[0] + uc*ka[0] + vc*kb[0],
                        faceA_center[1] + uc*ka[1] + vc*kb[1],
                        faceA_center[2] + uc*ka[2] + vc*kb[2]);
              log_trace("    d2 = %f %f %f\n",
                        faceA_center[0] + ud*ka[0] + vd*kb[0],
                        faceA_center[1] + ud*ka[1] + vd*kb[1],
                        faceA_center[2] + ud*ka[2] + vd*kb[2]);
            }
          }

          if (dbg_pair) {
            log_trace("    f = %f\n", f);
            log_trace("    area += %f\n", f*area_kab);
          }
          area += f*area_kab;

        } // End of loop on edges of current face B

      } // End of loop on edges of current face A

      a_to_b_volume[ifaceB] = signAB * area;
      if (0) {//dbg_enabled) {
        log_trace("faceA %d ("PDM_FMT_G_NUM") faceB %d ("PDM_FMT_G_NUM"), volume = %20.16f (%3.3f%)\n",
                  faceA_id, face_ln_to_gn[0][faceA_id],
                  faceB_id, face_ln_to_gn[1][faceB_id],
                  a_to_b_volume[ifaceB],
                  100*a_to_b_volume[ifaceB]*imag*2);
      }

    }  // End of loop on faces B

  } // End of loop on faces A


  if (dbg_enabled) {
    // Crude check
    double l_total_area_AB = 0;
    for (int i = 0; i < a_to_b_idx[n_face[0]]; i++) {
      l_total_area_AB += a_to_b_volume[i];
    }

    double l_total_area_A = 0;
    int idx = 0;
    for (int faceA_id = 0; faceA_id < n_face[0]; faceA_id++) {
      double faceA_normal[3];
      double faceA_center[3];
      if (face_vtx[0] == NULL) {
        _polygon_geom_properties(face_vtx_idx[0][faceA_id+1] - face_vtx_idx[0][faceA_id],
                                 face_edge[0] + face_vtx_idx[0][faceA_id],
                                 edge_vtx[0],
                                 vtx_coord[0],
                                 faceA_normal,
                                 faceA_center);
      }
      else {
        int is_degenerate = 0;
        PDM_geom_elem_polygon_properties(1,
                                         face_vtx_idx[0] + faceA_id,
                                         face_vtx[0],
                                         vtx_coord[0],
                                         faceA_normal,
                                         faceA_center,
                                         NULL,
                                         &is_degenerate);

        if (is_degenerate) {
          continue;
        }
      }
      double area = PDM_MODULE(faceA_normal);
      if (face_vtx[0] == NULL) {
        area *= 0.5;
      }

      if (1) {//faceA_ln_to_gn[faceA_id] == 2385) {
        double sum = 0;
        for (int j = a_to_b_idx[faceA_id]; j < a_to_b_idx[faceA_id+1]; j++) {
        // for (int j = 0; j < a_to_b_n[faceA_id]; j++) {
          log_trace(PDM_FMT_G_NUM"-"PDM_FMT_G_NUM" : %20.16f\n",
                    face_ln_to_gn[0][faceA_id], face_ln_to_gn[1][a_to_b[idx]], a_to_b_volume[idx]);
          sum += a_to_b_volume[idx];
          idx++;
        }
        log_trace(PDM_FMT_G_NUM" : sum = %20.16f / %20.16f (%f%%)\n",
                  face_ln_to_gn[0][faceA_id], sum,area, 100*sum/area);
      }

      l_total_area_A += area;
    }

    double g_total_area_AB;
    PDM_MPI_Allreduce(&l_total_area_AB, &g_total_area_AB, 1,
                      PDM_MPI_DOUBLE, PDM_MPI_SUM, mi->comm);

    double g_total_area_A;
    PDM_MPI_Allreduce(&l_total_area_A, &g_total_area_A, 1,
                      PDM_MPI_DOUBLE, PDM_MPI_SUM, mi->comm);

    log_trace("total area of A inter B : local = %20.16f, global = %20.16f (%3.3f%%)\n",
              l_total_area_AB, g_total_area_AB,
              100*g_total_area_AB / g_total_area_A);

    // cas plan, translation (0.5,0.5,0) + rotation PI/5
    double exact = 0.0875401518835469;
    log_trace("error : absolute = %e, relative = %e\n",
              PDM_ABS(g_total_area_AB - exact),
              PDM_ABS(g_total_area_AB - exact)/exact);

  }
  free(faceB_normals);
  for (int i = 0; i < 2; i++) {
    if (mi->mesh[i] == NULL) {
      free(face_vtx_idx[i]);
      free(face_vtx    [i]);
    }
  }

  /* Build part_to_part object */
  _build_ptp(mi,
             a_to_b_idx,
             a_to_b,
             a_to_b_volume);

  free(a_to_b_volume);
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_mesh_intersection_t*
PDM_mesh_intersection_create
(
 const PDM_mesh_intersection_kind_t intersection_kind,
 const int                          dim_mesh_a,
 const int                          dim_mesh_b,
 const double                       project_coeff,
       PDM_MPI_Comm                 comm,
 const PDM_ownership_t              owner
)
{
  PDM_mesh_intersection_t *mi = (PDM_mesh_intersection_t *) malloc(sizeof(PDM_mesh_intersection_t));

  mi->comm = comm;
  mi->intersect_kind = intersection_kind;
  mi->n_part_mesh[0] = 0;
  mi->n_part_mesh[1] = 0;
  mi->dim_mesh[0]    = dim_mesh_a;
  mi->dim_mesh[1]    = dim_mesh_b;
  mi->project_coef   = project_coeff;

  mi->bbox_tolerance = 1e-6;

  mi->mesh      [0] = NULL;
  mi->mesh      [1] = NULL;
  mi->mesh_nodal[0] = NULL;
  mi->mesh_nodal[1] = NULL;

  mi->extrp_mesh[0] = NULL;
  mi->extrp_mesh[1] = NULL;

  /* Initialize results */
  mi->owner = owner;
  mi->tag_extrp_mesh      = 0;
  mi->tag_elt_a_elt_b_get = 0;
  mi->tag_box_a_box_b_get = 0;
  mi->tag_elt_b_elt_a_get = 0;
  mi->elt_a_elt_b_idx     = NULL;
  mi->elt_a_elt_b         = NULL;
  mi->box_a_box_b_idx     = NULL;
  mi->box_a_box_b         = NULL;
  mi->elt_a_elt_b_volume  = NULL;
  mi->elt_b_elt_a_volume  = NULL;
  mi->ptp                 = NULL;
  mi->ptp_ownership       = PDM_OWNERSHIP_KEEP;

  mi->tag_elt_volume_get[0] = 0;
  mi->tag_elt_volume_get[1] = 0;
  mi->elt_volume[0] = NULL;
  mi->elt_volume[1] = NULL;


  return mi;
}

/**
 *
 * \brief Set global data of a mesh
 *
 * \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
 * \param [in]   i_mesh         Mesh identifier
 * \param [in]   n_part         Number of partitions
 *
 */

void
PDM_mesh_intersection_n_part_set
(
  PDM_mesh_intersection_t *mi,
  const int                i_mesh,
  const int                n_part
)
{
  assert(mi->mesh_nodal[i_mesh] == NULL);

  mi->n_part_mesh[i_mesh] = n_part;

  mi->mesh[i_mesh] = PDM_part_mesh_create(n_part, mi->comm);

  mi->elt_volume[i_mesh] = malloc(sizeof(double *) * mi->n_part_mesh[i_mesh]);
  for (int i = 0; i < mi->n_part_mesh[i_mesh]; i++) {
    mi->elt_volume[i_mesh][i] = NULL;
  }
}


void
PDM_mesh_intersection_compute
(
  PDM_mesh_intersection_t  *mi
)
{
  int i_rank;
  PDM_MPI_Comm_rank(mi->comm, &i_rank);

  /*
   * Compute extents of mesh_a and mesh_b
   */
  double **extents_mesh[2] = {NULL, NULL};
  double global_extents[6] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
  double mesh_global_extents[2][6] = {{ HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                                       -HUGE_VAL, -HUGE_VAL, -HUGE_VAL},
                                      {HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                                       -HUGE_VAL, -HUGE_VAL, -HUGE_VAL}};
  double g_mesh_global_extents[2][6];
  double g_global_extents        [6];
  for (int imesh = 0; imesh < 2; imesh++) {
    if (mi->mesh_nodal[imesh] == NULL && mi->mesh[imesh] != NULL) {
      _compute_part_mesh_extents(mi->mesh[imesh],
                                 mi->dim_mesh[imesh],
                                 mi->bbox_tolerance,
                                 mesh_global_extents[imesh],
                                 &extents_mesh[imesh]);
      // printf("extents_mesh [%i] = (%12.5e/%12.5e/%12.5e) | (%12.5e/%12.5e/%12.5e) \n", imesh,
      //        mesh_global_extents[imesh][0], mesh_global_extents[imesh][1], mesh_global_extents[imesh][2],
      //        mesh_global_extents[imesh][3], mesh_global_extents[imesh][4], mesh_global_extents[imesh][5]);
    }
    else if (mi->mesh_nodal[imesh] != NULL) {
      _compute_mesh_nodal_extents(mi->mesh_nodal[imesh],
                                  mi->dim_mesh[imesh],
                                  mi->bbox_tolerance,
                                  mesh_global_extents[imesh],
                                  &extents_mesh[imesh]);
    }
    else {
      extents_mesh[imesh] = malloc(sizeof(double *) * mi->n_part_mesh[imesh]);
      for (int i = 0; i < mi->n_part_mesh[imesh]; i++) {
        extents_mesh[imesh][i] = malloc(sizeof(double) * 0);
      }
    }
  }

  /*
   * Compute vertex normals if necessary
   */
  if (mi->dim_mesh[0] == 2 && mi->dim_mesh[1] == 2) {
    //...
  }



  /*
   * Global extents exchange
   */
  const int dim = 3;
  for(int i_mesh = 0; i_mesh < 2; ++i_mesh) {
    PDM_MPI_Allreduce(mesh_global_extents[i_mesh], g_mesh_global_extents[i_mesh], dim,
                      PDM_MPI_DOUBLE, PDM_MPI_MIN, mi->comm);
    PDM_MPI_Allreduce(mesh_global_extents[i_mesh]+dim, g_mesh_global_extents[i_mesh]+dim, dim,
                      PDM_MPI_DOUBLE, PDM_MPI_MAX, mi->comm);
  }

  /* Union or intersection of global extents */
  for(int i_mesh = 0; i_mesh < 2; ++i_mesh) {
    for (int k = 0; k < 3; k++) {
      // Union
      // global_extents[k]     = PDM_MIN(mesh_global_extents[i_mesh][k  ], global_extents[k  ]);
      // global_extents[3 + k] = PDM_MAX(mesh_global_extents[i_mesh][3+k], global_extents[3+k]);
      // Intersection
      global_extents[k]     = PDM_MAX(g_mesh_global_extents[i_mesh][k  ], global_extents[k  ]);
      global_extents[3 + k] = PDM_MIN(g_mesh_global_extents[i_mesh][3+k], global_extents[3+k]);
    }
  }
  for(int i = 0; i < 6; ++i) {
    g_global_extents[i] = global_extents[i];
  }
  double max_range = -HUGE_VAL;
  double min_range =  HUGE_VAL;

  for (int k = 0; k < dim; k++) {
    max_range = PDM_MAX(max_range, (g_global_extents[dim+k] - g_global_extents[k]));
    min_range = PDM_MIN(min_range, (g_global_extents[dim+k] - g_global_extents[k]));
  }

  for (int k = 0; k < dim; k++) {
    g_global_extents[k]     += -max_range * 1.1e-3; // On casse la symetrie !
    g_global_extents[dim+k] +=  max_range * 1e-3;
  }


  /*
   * Extraction - In option ?
   */
  int           n_mesh = 2;
  int           n_part                    [n_mesh];
  int          *n_extract_elmt            [n_mesh];
  double      **extract_box_extents       [n_mesh];
  int         **extract_elmt_init_location[n_mesh];
  PDM_g_num_t **extract_elmt_ln_to_gn     [n_mesh];

  n_part[0] = mi->n_part_mesh[0];
  n_part[1] = mi->n_part_mesh[1];
  for (int imesh = 0; imesh < 2; imesh++) {
    if (mi->mesh_nodal[imesh] == NULL && mi->mesh[imesh] != NULL) {
      _select_elements_by_global_bbox(mi->mesh[imesh],
                                      mi->dim_mesh[imesh],
                                      extents_mesh[imesh],
                                      g_mesh_global_extents[(imesh+1)%2],
                                      &n_extract_elmt[imesh],
                                      &extract_box_extents[imesh],
                                      &extract_elmt_init_location[imesh],
                                      &extract_elmt_ln_to_gn[imesh]);
    }
    else if (mi->mesh_nodal[imesh] != NULL) {
      _select_elements_by_global_bbox_nodal(mi,
                                            imesh,
                                            extents_mesh[imesh],
                                            g_mesh_global_extents[(imesh+1)%2],
                                            &n_extract_elmt[imesh],
                                            &extract_box_extents[imesh],
                                            &extract_elmt_init_location[imesh],
                                            &extract_elmt_ln_to_gn[imesh]);
    }
    else {
      n_extract_elmt            [imesh] = malloc(sizeof(int          ) * n_part[imesh]);
      extract_box_extents       [imesh] = malloc(sizeof(double      *) * n_part[imesh]);
      extract_elmt_init_location[imesh] = malloc(sizeof(int         *) * n_part[imesh]);
      extract_elmt_ln_to_gn     [imesh] = malloc(sizeof(PDM_g_num_t *) * n_part[imesh]);
      for (int ipart = 0; ipart < n_part[imesh]; ipart++) {
        n_extract_elmt[imesh][ipart] = 0;
        extract_box_extents       [imesh][ipart] = malloc(sizeof(double     ) * 0);
        extract_elmt_init_location[imesh][ipart] = malloc(sizeof(int        ) * 0);
        extract_elmt_ln_to_gn     [imesh][ipart] = malloc(sizeof(PDM_g_num_t) * 0);
      }
    }

    if (0  == 1) {
      for(int i_part = 0; i_part < n_part[imesh]; ++i_part) {
        log_trace("n_extract_elmt[%d][%d] = %d\n", imesh, i_part, n_extract_elmt[imesh][i_part]);
        // for (int i = 0; i < n_extract_elmt[imesh][i_part]; i++) {
        //   log_trace(" init_loc[%d] : (%d, %d, %d)\n", i,
        //             extract_elmt_init_location[imesh][i_part][3*i  ],
        //             extract_elmt_init_location[imesh][i_part][3*i+1],
        //             extract_elmt_init_location[imesh][i_part][3*i+2]);
        // }
        char filename[99];
        sprintf(filename, "select_boxes_mesh%d_part%d_rank%d.vtk", imesh, i_part, i_rank);
        PDM_vtk_write_boxes(filename,
                            n_extract_elmt[imesh][i_part],
                            extract_box_extents[imesh][i_part],
                            extract_elmt_ln_to_gn[imesh][i_part]);
      }
    }
  }

  for(int i_part = 0; i_part < n_part[0]; ++i_part) {
    free(extents_mesh[0][i_part]);
  }
  for(int i_part = 0; i_part < n_part[1]; ++i_part) {
    free(extents_mesh[1][i_part]);
  }
  free(extents_mesh[0]);
  free(extents_mesh[1]);

  // Attention le dbtree fait  le init_location sauf que la il faut le forcer !!!!!!
  PDM_dbbtree_t *dbbtree_mesh_a = PDM_dbbtree_create (mi->comm, dim, g_global_extents);

  PDM_box_set_t *boxes_mesh[2] = {NULL, NULL};
  boxes_mesh[0] = PDM_dbbtree_boxes_set_with_init_location(dbbtree_mesh_a,
                                                           n_part[0],
                                                           n_extract_elmt            [0],
                                   (const int         **)  extract_elmt_init_location[0],
                                   (const double      **)  extract_box_extents       [0],
                                   (const PDM_g_num_t **)  extract_elmt_ln_to_gn     [0]);

  /*
   * Intersect with B
   */
  int *box_a_to_box_b_idx = NULL;
  int *box_a_to_box_b     = NULL;
  boxes_mesh[1] = PDM_dbbtree_intersect_boxes_with_init_location_set(dbbtree_mesh_a,
                                                                     n_part[1],
                                                                     n_extract_elmt            [1],
                                              (const int         **) extract_elmt_init_location[1],
                                              (const double      **) extract_box_extents       [1],
                                              (const PDM_g_num_t **) extract_elmt_ln_to_gn     [1],
                                                                     &box_a_to_box_b_idx,
                                                                     &box_a_to_box_b);

  /* Free extraction */
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    for(int i_part = 0; i_part < n_part[i_mesh]; ++i_part) {
      free(extract_elmt_init_location[i_mesh][i_part]);
      free(extract_box_extents       [i_mesh][i_part]);
      free(extract_elmt_ln_to_gn     [i_mesh][i_part]);
    }
    free(n_extract_elmt            [i_mesh]);
    free(extract_elmt_init_location[i_mesh]);
    free(extract_box_extents       [i_mesh]);
    free(extract_elmt_ln_to_gn     [i_mesh]);
  }


  /*
   *  Redistrib all boxes (inplace) like overlay before extracting mesh
   */
  int *redistribute_box_a_to_box_b_idx = NULL;
  int *redistribute_box_a_to_box_b     = NULL;
  _redistrib_boxes(mi->comm,
                   boxes_mesh[0],
                   boxes_mesh[1],
                   box_a_to_box_b_idx,
                   box_a_to_box_b,
                   &redistribute_box_a_to_box_b_idx,
                   &redistribute_box_a_to_box_b);
  
  free(box_a_to_box_b_idx);
  free(box_a_to_box_b);


  /*
   * Extract part
   */
  // PDM_extract_part_t* extrp_mesh_a = _create_extract_part(mi->mesh[0],
  //                                                         mi->dim_mesh[0],
  //                                                         boxes_mesh[0]);
  // PDM_extract_part_t* extrp_mesh_b = _create_extract_part(mi->mesh[1],
  //                                                         mi->dim_mesh[1],
  //                                                         boxes_mesh[1]);
  for (int imesh = 0; imesh < 2; imesh++) {
    if (mi->mesh_nodal[imesh] == NULL && mi->mesh[imesh] != NULL) {
      mi->extrp_mesh[imesh] = _create_extract_part(mi->mesh[imesh],
                                                   mi->dim_mesh[imesh],
                                                   mi->intersect_kind,
                                                   boxes_mesh[imesh]);
    }
    else {
      mi->extrp_mesh[imesh] = _create_extract_part_nodal(mi,
                                                         imesh,
                                                         boxes_mesh[imesh]);
    }
  }
  PDM_extract_part_t* extrp_mesh_a = mi->extrp_mesh[0];
  PDM_extract_part_t* extrp_mesh_b = mi->extrp_mesh[1];

  PDM_dbbtree_free (dbbtree_mesh_a);
  // PDM_box_set_destroy (&boxes_mesh[0]);
  // PDM_box_set_destroy (&boxes_mesh[1]);

  /*
   * Geometry begin here ...
   */

  /* Visu */
  if (0 == 1) {
    for (int imesh = 0; imesh < 2; imesh++) {
      if (mi->mesh_nodal[imesh] != NULL) {

        PDM_part_mesh_nodal_t *extract_pmn = PDM_part_mesh_nodal_create(mi->dim_mesh[imesh],
                                                                        1,
                                                                        mi->comm);

        double *_vtx_coord = NULL;
        int _n_vtx = PDM_extract_part_vtx_coord_get(mi->extrp_mesh[imesh],
                                                    0,
                                                    &_vtx_coord,
                                                    PDM_OWNERSHIP_KEEP);
        PDM_g_num_t *_vtx_ln_to_gn = NULL;
        PDM_extract_part_parent_ln_to_gn_get(mi->extrp_mesh[imesh],
                                             0,
                                             PDM_MESH_ENTITY_VTX,
                                             &_vtx_ln_to_gn,
                                             PDM_OWNERSHIP_KEEP);

        PDM_part_mesh_nodal_coord_set(extract_pmn,
                                      0,
                                      _n_vtx,
                                      _vtx_coord,
                                      _vtx_ln_to_gn,
                                      PDM_OWNERSHIP_USER);

        PDM_part_mesh_nodal_elmts_t *extract_pmne = NULL;
        PDM_extract_part_part_mesh_nodal_get(mi->extrp_mesh[imesh],
                                             &extract_pmne,
                                             PDM_OWNERSHIP_USER);

        PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(extract_pmn, extract_pmne);


        char pattern[99];
        sprintf(pattern, "extrp_mesh%d", imesh);

        PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_MAX;
        switch (mi->dim_mesh[imesh]) {
        case 1:
          geom_kind = PDM_GEOMETRY_KIND_RIDGE;
          break;
        case 2:
          geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
          break;
        case 3:
          geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
          break;
        default:
          PDM_error(__FILE__, __LINE__, 0, "invalid dimension %d\n", mi->dim_mesh[imesh]);
        }

        PDM_part_mesh_nodal_dump_vtk(extract_pmn,
                                     geom_kind,
                                     pattern);

        PDM_part_mesh_nodal_free(extract_pmn);
      }
    }
  }


  PDM_MPI_Barrier(mi->comm);

  if (mi->intersect_kind == PDM_MESH_INTERSECTION_KIND_PREPROCESS) {
    mi->box_a_box_b_idx   = redistribute_box_a_to_box_b_idx;
    mi->box_a_box_b       = redistribute_box_a_to_box_b;
    PDM_box_set_destroy (&boxes_mesh[0]);
    PDM_box_set_destroy (&boxes_mesh[1]);

  }

  else {

    // if (mi->mesh_nodal[0] != NULL ||
    //     mi->mesh_nodal[1] != NULL) {
    //   PDM_error(__FILE__, __LINE__, 0, "Nodal version not implemented yet\n");
    // }

    if(mi->dim_mesh[0] == 3 && mi->dim_mesh[1] == 3) {
        _mesh_intersection_vol_vol(mi,
                                   redistribute_box_a_to_box_b_idx,
                                   redistribute_box_a_to_box_b);
    } else if(mi->dim_mesh[0] == 3 && mi->dim_mesh[1] == 2) {
      // On suppose que l'utilisateur met A = Vol et B = Surf
      _mesh_intersection_vol_surf(mi,
                                  extrp_mesh_a,
                                  extrp_mesh_b,
                                  redistribute_box_a_to_box_b_idx,
                                  redistribute_box_a_to_box_b);
    } else if(mi->dim_mesh[0] == 2 && mi->dim_mesh[1] == 2) {
      _mesh_intersection_surf_surf(mi,
                                   redistribute_box_a_to_box_b_idx,
                                   redistribute_box_a_to_box_b);
    } else if(mi->dim_mesh[0] == 2 && mi->dim_mesh[1] == 1) {
      // On suppose que l'utilisateur met A = Vol et B = Surf
      _mesh_intersection_surf_line(mi,
                                   extrp_mesh_a,
                                   extrp_mesh_b,
                                   redistribute_box_a_to_box_b_idx,
                                   redistribute_box_a_to_box_b);
    } else if(mi->dim_mesh[0] == 3 && mi->dim_mesh[1] == 1) {
      // On suppose que l'utilisateur met A = Vol et B = Surf
      _mesh_intersection_vol_line(mi,
                                  extrp_mesh_a,
                                  extrp_mesh_b,
                                  redistribute_box_a_to_box_b_idx,
                                  redistribute_box_a_to_box_b);
    } else {
      PDM_error(__FILE__, __LINE__, 0,
                "PDM_mesh_intersection_compute error : Cannot handle meshA with dim = %i and meshB = %i \n", mi->dim_mesh[0], mi->dim_mesh[1]);
    }
    PDM_box_set_destroy (&boxes_mesh[0]);
    PDM_box_set_destroy (&boxes_mesh[1]);

    free(redistribute_box_a_to_box_b_idx);
    free(redistribute_box_a_to_box_b    );

    for (int imesh = 0; imesh < 2; imesh++) {
      if (mi->mesh_nodal[imesh] == NULL && mi->mesh[imesh] == NULL) {
        PDM_part_mesh_nodal_elmts_free(mi->extrp_mesh[imesh]->pmne);
      }
    }

    PDM_extract_part_free(extrp_mesh_a);
    PDM_extract_part_free(extrp_mesh_b);

    mi->extrp_mesh[0] = NULL; 
    mi->extrp_mesh[1] = NULL;
  }

}

/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   mi             Pointer to \ref PDM_mesh_intersection object
 * \param [in]   i_mesh         Mesh identifier
 * \param [in]   mesh           Pointer to \ref PDM_part_mesh_nodal object
 *
 */

void
PDM_mesh_intersection_mesh_nodal_set
(
 PDM_mesh_intersection_t  *mi,
 int                       i_mesh,
 PDM_part_mesh_nodal_t    *mesh
 )
{
  mi->mesh_nodal[i_mesh] = mesh;

  mi->n_part_mesh[i_mesh] = PDM_part_mesh_nodal_n_part_get(mesh);

  mi->elt_volume[i_mesh] = malloc(sizeof(double *) * mi->n_part_mesh[i_mesh]);
  for (int i = 0; i < mi->n_part_mesh[i_mesh]; i++) {
    mi->elt_volume[i_mesh][i] = NULL;
  }
}


void
PDM_mesh_intersection_part_set
(
  PDM_mesh_intersection_t  *mi,
  int                       i_mesh,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
)
{
  assert(i_part < mi->n_part_mesh[i_mesh]);

  PDM_part_mesh_t* mesh = mi->mesh[i_mesh];

  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_CELL, n_cell);
  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_FACE, n_face);
  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_EDGE, n_edge);
  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_VTX,  n_vtx );

  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, cell_face, cell_face_idx, PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE, face_edge, face_edge_idx, PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , face_vtx , face_vtx_idx , PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , edge_vtx , NULL         , PDM_OWNERSHIP_USER);

  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_CELL, cell_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_FACE, face_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_EDGE, edge_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_VTX,  vtx_ln_to_gn , PDM_OWNERSHIP_USER);

  PDM_part_mesh_vtx_coord_set(mesh, i_part, vtx_coord, PDM_OWNERSHIP_USER);
}

/* vol_vol */
void
PDM_mesh_intersection_tetraisation_pt_set
(
 PDM_mesh_intersection_t* mi,
 int     tetraisation_pt_type,
 double *tetraisation_pt_coord
)
{
  mi->tetraisation_pt_type = tetraisation_pt_type;
  // mi->tetraisation_pt_coord = NULL;
  if (tetraisation_pt_type == 2) {
    // mi->tetraisation_pt_coord = malloc(sizeof(double) * 3);
    memcpy(mi->tetraisation_pt_coord, tetraisation_pt_coord, sizeof(double) * 3);
  }
}

/* debug */
void
PDM_mesh_intersection_stat_get
(
 PDM_mesh_intersection_t* mi,
 double *local_vol_A_B,
 double *global_vol_A_B,
 double *global_vol_A

)
{
  *local_vol_A_B  = mi->local_vol_A_B;
  *global_vol_A_B = mi->global_vol_A_B;
  *global_vol_A   = mi->global_vol_A;
}

void
PDM_mesh_intersection_free
(
 PDM_mesh_intersection_t* mi
)
{

  PDM_part_mesh_free(mi->mesh[0]);
  PDM_part_mesh_free(mi->mesh[1]);


  if (mi->ptp != NULL && mi->ptp_ownership == PDM_OWNERSHIP_KEEP) {
    PDM_part_to_part_free(mi->ptp);
  }

  if (mi->elt_a_elt_b_idx != NULL) {
    if ((mi->owner == PDM_OWNERSHIP_KEEP ) ||
        (mi->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !mi->tag_elt_a_elt_b_get)) {
      for (int ipart = 0; ipart < mi->n_part_mesh[0]; ipart++) {
        if (mi->elt_a_elt_b_idx[ipart] != NULL) {
          free(mi->elt_a_elt_b_idx[ipart]);
        }
      }
    }
    free(mi->elt_a_elt_b_idx);
  }

  if (mi->elt_a_elt_b != NULL) {
    if ((mi->owner == PDM_OWNERSHIP_KEEP ) ||
        (mi->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !mi->tag_elt_a_elt_b_get)) {
      for (int ipart = 0; ipart < mi->n_part_mesh[0]; ipart++) {
        if (mi->elt_a_elt_b[ipart] != NULL) {
          free(mi->elt_a_elt_b[ipart]);
        }
      }
  }
    free(mi->elt_a_elt_b);
  }

  if (mi->box_a_box_b != NULL) {
    if ((mi->owner == PDM_OWNERSHIP_KEEP ) ||
        (mi->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !mi->tag_box_a_box_b_get)) {
      free(mi->box_a_box_b);
      mi->box_a_box_b = NULL;
    }
  }

  if (mi->box_a_box_b_idx != NULL) {
    if ((mi->owner == PDM_OWNERSHIP_KEEP ) ||
        (mi->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !mi->tag_box_a_box_b_get)) {
      free(mi->box_a_box_b_idx);
      mi->box_a_box_b = NULL;
    }
  }

  if ((mi->owner == PDM_OWNERSHIP_KEEP ) ||
      (mi->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !mi->tag_extrp_mesh)) {
    for (int imesh = 0; imesh < 2; imesh++) {
      if (mi->extrp_mesh[imesh] != NULL) {
        PDM_extract_part_free(mi->extrp_mesh[imesh]);
        mi->extrp_mesh[imesh] = NULL;
      }
    }
  }

  if (mi->elt_a_elt_b_volume != NULL) {
    if ((mi->owner == PDM_OWNERSHIP_KEEP ) ||
        (mi->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !mi->tag_elt_a_elt_b_get)) {
      for (int ipart = 0; ipart < mi->n_part_mesh[0]; ipart++) {
        if (mi->elt_a_elt_b_volume[ipart] != NULL) {
          free(mi->elt_a_elt_b_volume[ipart]);
        }
      }
    }
    free(mi->elt_a_elt_b_volume);
  }

  if (mi->elt_b_elt_a_volume != NULL) {
    if ((mi->owner == PDM_OWNERSHIP_KEEP ) ||
        (mi->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !mi->tag_elt_b_elt_a_get)) {
      for (int ipart = 0; ipart < mi->n_part_mesh[1]; ipart++) {
        if (mi->elt_b_elt_a_volume[ipart] != NULL) {
          free(mi->elt_b_elt_a_volume[ipart]);
        }
      }
    }
    free(mi->elt_b_elt_a_volume);
  }

  for (int imesh = 0; imesh < 2; imesh++) {
    if (mi->elt_volume[imesh] != NULL) {
      if ((mi->owner == PDM_OWNERSHIP_KEEP ) ||
          (mi->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !mi->tag_elt_volume_get[imesh])) {
        for (int ipart = 0; ipart < mi->n_part_mesh[imesh]; ipart++) {
          if (mi->elt_volume[imesh][ipart] != NULL) {
            free(mi->elt_volume[imesh][ipart]);
          }
        }
      }
      free(mi->elt_volume[imesh]);
    }
  }

  free(mi);
}

/**
 * \brief Get part_to_part object to exchange data between the intersected meshes
 *
 * \param [in ] mi         Pointer to \ref PDM_mesh_intersection_t object
 * \param [out] ptp        Pointer to \ref PDM_part_to_part_t object
 * \param [in ] ownership  Ownership for ptp
 *
 */

void
PDM_mesh_intersection_part_to_part_get
(
 PDM_mesh_intersection_t  *mi,
 PDM_part_to_part_t      **ptp,
 PDM_ownership_t           ownership
 )
{
  *ptp = mi->ptp;
  mi->ptp_ownership = ownership;
}



void
PDM_mesh_intersection_result_from_a_get
(
       PDM_mesh_intersection_t  *mi,
 const int                       ipart,
       int                     **elt_a_elt_b_idx,
       PDM_g_num_t             **elt_a_elt_b,
       double                  **elt_a_elt_b_volume
 )
{
  assert(ipart < mi->n_part_mesh[0]);

  *elt_a_elt_b_idx    = mi->elt_a_elt_b_idx   [ipart];
  *elt_a_elt_b        = mi->elt_a_elt_b       [ipart];
  *elt_a_elt_b_volume = mi->elt_a_elt_b_volume[ipart];

  mi->tag_elt_a_elt_b_get = 1;
}


void
PDM_mesh_intersection_result_from_b_get
(
       PDM_mesh_intersection_t  *mi,
 const int                       ipart,
       double                  **elt_b_elt_a_volume
 )
{
  assert(ipart < mi->n_part_mesh[1]);

  *elt_b_elt_a_volume = mi->elt_b_elt_a_volume[ipart];

  mi->tag_elt_b_elt_a_get = 1;
}


void
PDM_mesh_intersection_elt_volume_get
(
       PDM_mesh_intersection_t  *mi,
 const int                       imesh,
 const int                       ipart,
       double                  **elt_volume
 )
{
  assert(imesh >= 0 && imesh < 2);
  assert(ipart < mi->n_part_mesh[imesh]);

  if (mi->elt_volume[imesh][ipart] == NULL) {
    /* Compute elt volumes */
    if (mi->mesh_nodal[imesh] == NULL && mi->mesh[imesh] != NULL) {

      int n_vtx = PDM_part_mesh_n_entity_get(mi->mesh[imesh],
                                             ipart,
                                             PDM_MESH_ENTITY_VTX);
      double *vtx_coord = NULL;
      PDM_part_mesh_vtx_coord_get(mi->mesh[imesh],
                                  ipart,
                                  &vtx_coord,
                                  PDM_OWNERSHIP_BAD_VALUE);

      int n_face = PDM_part_mesh_n_entity_get(mi->mesh[imesh],
                                              ipart,
                                              PDM_MESH_ENTITY_FACE);
      int *face_vtx_idx = NULL;
      int *face_vtx     = NULL;
      PDM_part_mesh_connectivity_get(mi->mesh[imesh],
                                     ipart,
                                     PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                     &face_vtx,
                                     &face_vtx_idx,
                                     PDM_OWNERSHIP_BAD_VALUE);

      int n_edge = PDM_part_mesh_n_entity_get(mi->mesh[imesh],
                                              ipart,
                                              PDM_MESH_ENTITY_EDGE);
      int *edge_vtx_idx = NULL;
      int *edge_vtx     = NULL;
      PDM_part_mesh_connectivity_get(mi->mesh[imesh],
                                     ipart,
                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                     &edge_vtx,
                                     &edge_vtx_idx,
                                     PDM_OWNERSHIP_BAD_VALUE);

      int *_face_vtx = face_vtx;
      if (mi->dim_mesh[imesh] > 1 && face_vtx == NULL) {
        int *face_edge = NULL;
        PDM_part_mesh_connectivity_get(mi->mesh[imesh],
                                       ipart,
                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                       &face_edge,
                                       &face_vtx_idx,
                                       PDM_OWNERSHIP_BAD_VALUE);

        PDM_compute_face_vtx_from_face_and_edge(n_face,
                                                face_vtx_idx,
                                                face_edge,
                                                edge_vtx,
                                                &_face_vtx);
      }

      switch (mi->dim_mesh[imesh]) {
      case 1: {
       mi->elt_volume[imesh][ipart] = malloc(sizeof(double) * n_edge);
       double *center = malloc(sizeof(double) * n_face * 3);
       PDM_geom_elem_edges_properties(n_edge,
                                      edge_vtx,
                                      vtx_coord,
                                      mi->elt_volume[imesh][ipart],
                                      center,
                                      NULL,
                                      NULL);
       free(center);
        break;
      }
      case 2: {
        mi->elt_volume[imesh][ipart] = malloc(sizeof(double) * n_face * 3);
        double *center = malloc(sizeof(double) * n_face * 3);
        PDM_geom_elem_polygon_properties(n_face,
                                         face_vtx_idx,
                                         _face_vtx,
                                         vtx_coord,
                                         mi->elt_volume[imesh][ipart],
                                         center,
                                         NULL,
                                         NULL);
        free(center);
        for (int i = 0; i < n_face; i++) {
          mi->elt_volume[imesh][ipart][i] = PDM_MODULE(mi->elt_volume[imesh][ipart] + 3*i);
        }
        mi->elt_volume[imesh][ipart] = realloc(mi->elt_volume[imesh][ipart],
                                               sizeof(double) * n_face);
        break;
      }
      case 3: {
        int n_cell = PDM_part_mesh_n_entity_get(mi->mesh[imesh],
                                                ipart,
                                                PDM_MESH_ENTITY_CELL);
        int *cell_face     = NULL;
        int *cell_face_idx = NULL;
        PDM_part_mesh_connectivity_get(mi->mesh[imesh],
                                       ipart,
                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                       &cell_face,
                                       &cell_face_idx,
                                       PDM_OWNERSHIP_BAD_VALUE);




        mi->elt_volume[imesh][ipart] = malloc(sizeof(double) * n_cell);
        double *center = malloc(sizeof(double) * n_cell * 3);
        PDM_geom_elem_polyhedra_properties_triangulated(1,
                                                        n_cell,
                                                        n_face,
                                                        face_vtx_idx,
                                                        _face_vtx,
                                                        cell_face_idx,
                                                        cell_face,
                                                        n_vtx,
                                                        vtx_coord,
                                                        mi->elt_volume[imesh][ipart],
                                                        center,
                                                        NULL,
                                                        NULL);
        free(center);
        break;
      }
      default:
        PDM_error(__FILE__, __LINE__, 0, "Invalid dim_mesh %d for mesh %d\n", mi->dim_mesh[imesh], imesh);
      }



      if (mi->dim_mesh[imesh] > 1 && face_vtx == NULL) {
        free(_face_vtx);
      }
    }
    else if (mi->mesh_nodal[imesh] != NULL) {
      // TODO...
      PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for nodal mesh\n");
    }
  }

  *elt_volume = mi->elt_volume[imesh][ipart];

  mi->tag_elt_volume_get[imesh] = 1;
}



/**
 *
 * \brief Set the tolerance for bounding boxes
 *
 * \param [in]   mi              Pointer to \ref PDM_mesh_intersection object
 * \param [in]   tol             Tolerance
 *
 */
void
PDM_mesh_intersection_tolerance_set
(
       PDM_mesh_intersection_t *mi,
 const double                   tol
)
{

  mi->bbox_tolerance = tol;
}


/**
 * \brief Get preprocessing results 
 *
 * \param [in ] mi                 Pointer to \ref PDM_mesh_intersection_t object
 * \param [out] elt_a_elt_b_idx    Index of list of intersected B element candidate for each A element
 *                                 in the extr_mesh distribution 
 * \param [out] elt_a_elt_b        List of intersected B element candidate for each A element in the 
 *                                 extr_mesh distribution 
 * \param [out]                    Redistributed mesh A with only A element candidate  
 * \param [out]                    Redistributed mesh B with only B element candidate  
 *
 */

void
PDM_mesh_intersection_preprocessing_get
(
       PDM_mesh_intersection_t  *mi,
       int                     **box_a_box_b_idx,
       int                     **box_a_box_b,
       PDM_extract_part_t      **extr_mesh_a,
       PDM_extract_part_t      **extr_mesh_b
)
{

  assert(mi->intersect_kind == PDM_MESH_INTERSECTION_KIND_PREPROCESS);
  *extr_mesh_a = mi->extrp_mesh[0];
  *extr_mesh_b = mi->extrp_mesh[1];
   mi->tag_extrp_mesh = 1;
  *box_a_box_b_idx = mi->box_a_box_b_idx;
  *box_a_box_b = mi->box_a_box_b;  
}


/**
 *
 * \brief Get mesh dimension
 *
 * \param [in] mi                 Pointer to \ref PDM_mesh_intersection_t object
 * \param [in] i_mesh             Mesh identifier
 *
 * \return Dimension of mesh \p i_mesh
 */

int
PDM_mesh_intersection_mesh_dimension_get
(
       PDM_mesh_intersection_t  *mi,
 const int                       i_mesh
)
{
  assert(i_mesh == 0 || i_mesh == 1);
  return mi->dim_mesh[i_mesh];
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
