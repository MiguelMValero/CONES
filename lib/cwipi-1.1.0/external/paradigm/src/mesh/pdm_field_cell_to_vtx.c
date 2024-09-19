
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_logging.h"
#include "pdm_distant_neighbor.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_order.h"
#include "pdm_binary_search.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_domain_interface_priv.h"
#include "pdm_field_cell_to_vtx_priv.h"
#include "pdm_field_cell_to_vtx.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/


/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

// static
// void
// triplet_to_quadruplet
// (
//   int   size,
//   int  *triplet,
//   int  *array,
//   int **quadruplet
// )
// {
//   int *_quadruplet = malloc(4 * size * sizeof(int));
//   for(int i = 0; i < size; ++i) {
//     _quadruplet[4*i  ] = triplet[3*i  ];
//     _quadruplet[4*i+1] = triplet[3*i+1];
//     _quadruplet[4*i+2] = triplet[3*i+2];
//     _quadruplet[4*i+3] = array  [i];
//   }


//   *quadruplet = _quadruplet;
// }

// static inline
// int
// _is_same_quadruplet
// (
// int iproc1, int ipart1, int ielt1, int iinterf1,
// int iproc2, int ipart2, int ielt2, int iinterf2
// )
// {
//   if(iproc1 == iproc2){
//     if(ipart1 == ipart2){
//       if(ielt1 == ielt2){
//         if(iinterf1 == iinterf2){
//           return 1;
//         }
//       }
//     }
//   }
//   return 0;
// }

static
void
_cell_center_3d
(
  int       n_part_in,
  int      *pn_cell,
  int     **pcell_face_idx,
  int     **pcell_face,
  int     **pface_edge_idx,
  int     **pface_edge,
  int     **pface_vtx_idx,
  int     **pface_vtx,
  int     **pedge_vtx,
  double  **pvtx_coord,
  double ***cell_center
)
{
  int from_edge = 0;
  int from_face = 0;
  for(int i_part = 0; i_part < n_part_in; ++i_part) {
    if(pface_edge    [i_part] != NULL) {
      from_edge = 1;
    }
    if(pface_vtx    [i_part] != NULL) {
      from_face = 1;
    }
    assert(pvtx_coord    [i_part] != NULL);
  }

  double** entity_center = malloc(n_part_in * sizeof(double * ));

  if(from_face == 1) {
    for(int i_part = 0; i_part < n_part_in; ++i_part) {
      entity_center[i_part] = (double *) malloc(3 * pn_cell[i_part] * sizeof(double));

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_vtx      = pface_vtx     [i_part];
      int    *_pface_vtx_idx  = pface_vtx_idx [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      // PDM_log_trace_array_int(extract_lnum[i_part], pn_cell[i_part], "extract_lnum ::");
      for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {

        entity_center[i_part][3*i_cell  ] = 0.;
        entity_center[i_part][3*i_cell+1] = 0.;
        entity_center[i_part][3*i_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[i_cell+1] - _pcell_face_idx[i_cell]);

        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double inv2 = 1./((double)  _pface_vtx_idx[i_face+1] - _pface_vtx_idx[i_face]);

          double fcx = 0;
          double fcy = 0;
          double fcz = 0;
          for(int idx_vtx = _pface_vtx_idx[i_face]; idx_vtx < _pface_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = _pface_vtx[idx_vtx]-1;
            fcx += _pvtx_coord[3*i_vtx  ];
            fcy += _pvtx_coord[3*i_vtx+1];
            fcz += _pvtx_coord[3*i_vtx+2];
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*i_cell  ] += fcx;
          entity_center[i_part][3*i_cell+1] += fcy;
          entity_center[i_part][3*i_cell+2] += fcz;
        }

        entity_center[i_part][3*i_cell  ] = entity_center[i_part][3*i_cell  ] * inv;
        entity_center[i_part][3*i_cell+1] = entity_center[i_part][3*i_cell+1] * inv;
        entity_center[i_part][3*i_cell+2] = entity_center[i_part][3*i_cell+2] * inv;
      } /* End cell */
    }
  } else if( from_edge == 1) {
    for(int i_part = 0; i_part < n_part_in; ++i_part) {
      entity_center[i_part] = (double *) malloc(3 * pn_cell[i_part] * sizeof(double));

      int    *_pcell_face     = pcell_face    [i_part];
      int    *_pcell_face_idx = pcell_face_idx[i_part];
      int    *_pface_edge     = pface_edge    [i_part];
      int    *_pface_edge_idx = pface_edge_idx[i_part];
      int    *_pedge_vtx      = pedge_vtx     [i_part];
      double *_pvtx_coord     = pvtx_coord    [i_part];

      for(int i_cell = 0; i_cell < pn_cell[i_part]; ++i_cell) {

        entity_center[i_part][3*i_cell  ] = 0.;
        entity_center[i_part][3*i_cell+1] = 0.;
        entity_center[i_part][3*i_cell+2] = 0.;

        double inv = 1./((double)  _pcell_face_idx[i_cell+1] - _pcell_face_idx[i_cell]);

        double fcx = 0;
        double fcy = 0;
        double fcz = 0;
        for(int idx_face = _pcell_face_idx[i_cell]; idx_face < _pcell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(_pcell_face[idx_face])-1;

          double inv2 = 1./((double)  _pface_edge_idx[i_face+1] - _pface_edge_idx[i_face]);

          for(int idx_edge = _pface_edge_idx[i_face]; idx_edge < _pface_edge_idx[i_face+1]; ++idx_edge) {
            int i_edge = PDM_ABS(_pface_edge[idx_edge])-1;
            int i_vtx1 = _pedge_vtx[2*i_edge  ] - 1;
            int i_vtx2 = _pedge_vtx[2*i_edge+1] - 1;
            fcx += 0.5 * (_pvtx_coord[3*i_vtx1  ] + _pvtx_coord[3*i_vtx2  ]);
            fcy += 0.5 * (_pvtx_coord[3*i_vtx1+1] + _pvtx_coord[3*i_vtx2+1]);
            fcz += 0.5 * (_pvtx_coord[3*i_vtx1+2] + _pvtx_coord[3*i_vtx2+2]);
          }
          fcx = fcx * inv2;
          fcy = fcy * inv2;
          fcz = fcz * inv2;

          entity_center[i_part][3*i_cell  ] += fcx;
          entity_center[i_part][3*i_cell+1] += fcy;
          entity_center[i_part][3*i_cell+2] += fcz;
        }

        entity_center[i_part][3*i_cell  ] = entity_center[i_part][3*i_cell  ] * inv;
        entity_center[i_part][3*i_cell+1] = entity_center[i_part][3*i_cell+1] * inv;
        entity_center[i_part][3*i_cell+2] = entity_center[i_part][3*i_cell+2] * inv;
      } /* End cell */
    }
  }

  *cell_center = entity_center;
}

static
void
_prepare_cell_center
(
  PDM_field_cell_to_vtx_t* fctv
)
{
  int     *pn_cell        = malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  int    **pcell_face_idx = malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  int    **pcell_face     = malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  int    **pface_edge_idx = malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  int    **pface_edge     = malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  int    **pface_vtx_idx  = malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  int    **pface_vtx      = malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  int    **pedge_vtx      = malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  double **pvtx_coord     = malloc(fctv->n_part_loc_all_domain * sizeof(double *));

  int shift_part = 0;
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
      pn_cell       [i_part+shift_part] = fctv->parts[i_domain][i_part].n_cell;
      pcell_face_idx[i_part+shift_part] = fctv->parts[i_domain][i_part].cell_face_idx;
      pcell_face    [i_part+shift_part] = fctv->parts[i_domain][i_part].cell_face;
      pface_edge_idx[i_part+shift_part] = fctv->parts[i_domain][i_part].face_edge_idx;
      pface_edge    [i_part+shift_part] = fctv->parts[i_domain][i_part].face_edge;
      pface_vtx_idx [i_part+shift_part] = fctv->parts[i_domain][i_part].face_vtx_idx;
      pface_vtx     [i_part+shift_part] = fctv->parts[i_domain][i_part].face_vtx;
      pedge_vtx     [i_part+shift_part] = fctv->parts[i_domain][i_part].edge_vtx;
      pvtx_coord    [i_part+shift_part] = fctv->parts[i_domain][i_part].vtx;
    }
    shift_part += fctv->n_part[i_domain];
  }

  assert(fctv->cell_center == NULL);
  _cell_center_3d(fctv->n_part_loc_all_domain,
                  pn_cell,
                  pcell_face_idx,
                  pcell_face,
                  pface_edge_idx,
                  pface_edge,
                  pface_vtx_idx,
                  pface_vtx,
                  pedge_vtx,
                  pvtx_coord,
                  &fctv->cell_center);

  if(0 == 1) {

    int i_rank;
    PDM_MPI_Comm_rank(fctv->comm, &i_rank);

    shift_part = 0;
    for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
        char filename[999];
        sprintf(filename, "vtx_coords_%i_%i.vtk", i_rank, i_part+shift_part);
        PDM_vtk_write_point_cloud(filename,
                                  fctv->parts[i_domain][i_part].n_vtx,
                                  fctv->parts[i_domain][i_part].vtx,
                                  NULL,
                                  NULL);

        sprintf(filename, "cell_center_%i_%i.vtk", i_rank, i_part+shift_part);
        PDM_vtk_write_point_cloud(filename,
                                  fctv->parts[i_domain][i_part].n_cell,
                                  fctv->cell_center[i_part+shift_part],
                                  NULL,
                                  NULL);
      }
      shift_part += fctv->n_part[i_domain];
    }
  }

  free(pn_cell);
  free(pcell_face_idx);
  free(pcell_face    );
  free(pface_edge_idx);
  free(pface_edge    );
  free(pface_vtx_idx );
  free(pface_vtx     );
  free(pedge_vtx     );
  free(pvtx_coord    );
}

static
void
_prepare_vtx_cell
(
 PDM_field_cell_to_vtx_t*  fctv
)
{
  /*
   * Create vtx_cell
   */
  int shift_part = 0;
  fctv->pvtx_cell_n   = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int *));
  fctv->pvtx_cell_idx = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int *));
  fctv->pvtx_cell     = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int *));
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {

      /* Compute cell_edge */
      int* cell_vtx_idx = NULL;
      int* cell_vtx     = NULL;
      PDM_combine_connectivity(fctv->parts[i_domain][i_part].n_cell,
                               fctv->parts[i_domain][i_part].cell_face_idx,
                               fctv->parts[i_domain][i_part].cell_face,
                               fctv->parts[i_domain][i_part].face_vtx_idx,
                               fctv->parts[i_domain][i_part].face_vtx,
                               &cell_vtx_idx,
                               &cell_vtx);

      PDM_connectivity_transpose(fctv->parts[i_domain][i_part].n_cell,
                                 fctv->parts[i_domain][i_part].n_vtx,
                                 cell_vtx_idx,
                                 cell_vtx,
                                 &fctv->pvtx_cell_idx[i_part+shift_part],
                                 &fctv->pvtx_cell    [i_part+shift_part]);

      free(cell_vtx_idx);
      free(cell_vtx);

      fctv->pvtx_cell_n[i_part+shift_part] = malloc(fctv->parts[i_domain][i_part].n_vtx * sizeof(int));
      for(int i_vtx = 0; i_vtx < fctv->parts[i_domain][i_part].n_vtx; ++i_vtx) {
        fctv->pvtx_cell_n[i_part+shift_part][i_vtx] = fctv->pvtx_cell_idx[i_part+shift_part][i_vtx+1] - fctv->pvtx_cell_idx[i_part+shift_part][i_vtx];
      }
    }
    shift_part += fctv->n_part[i_domain];
  }
}

static
void
_prepare_cell_center_nodal
(
  PDM_field_cell_to_vtx_t* fctv
)
{
  assert(fctv->cell_center == NULL);
  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
  int order = 1;

  int shift_part = 0;
  fctv->cell_center = malloc(fctv->n_part_loc_all_domain * sizeof(double * ));
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {

    PDM_part_mesh_nodal_elmts_t* pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(fctv->pmn[i_domain], geom_kind);

    int n_section    = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
    int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {

      // int n_vtx        = PDM_part_mesh_nodal_elmts_n_vtx_get      (pmne, i_part);

      int n_elmt = 0;
      for(int i_section = 0; i_section < n_section; ++i_section) {
        int n_elmt_in_section = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                                            sections_id[i_section],
                                                                            i_part);
        n_elmt     += n_elmt_in_section;
      }

      fctv->cell_center[i_part+shift_part] = malloc(3 * n_elmt * sizeof(double));
      double *_entity_center = fctv->cell_center[i_part+shift_part];
      double *pvtx_coord = fctv->parts[i_domain][i_part].vtx;

      n_elmt = 0;
      for(int i_section = 0; i_section < n_section; ++i_section) {
        PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, sections_id[i_section]);
        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
        int n_elmt_in_section = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                                            sections_id[i_section],
                                                                            i_part);
        int            *connec              = NULL;
        PDM_g_num_t    *numabs              = NULL;
        int            *parent_num          = NULL;
        PDM_g_num_t    *parent_entity_g_num = NULL;
        PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                  sections_id[i_section],
                                                  i_part,
                                                  &connec,
                                                  &numabs,
                                                  &parent_num,
                                                  &parent_entity_g_num,
                                                  PDM_OWNERSHIP_KEEP);

        double pond = 1./((double) n_vtx_per_elmt);
        for(int i_elmt = 0; i_elmt < n_elmt_in_section; ++i_elmt) {
          _entity_center[3*n_elmt  ] = 0.;
          _entity_center[3*n_elmt+1] = 0.;
          _entity_center[3*n_elmt+2] = 0.;
          int i_vtx = connec[i_elmt]-1;
          for(int k = 0; k < n_vtx_per_elmt; ++k) {
            _entity_center[3*n_elmt  ] += pvtx_coord[3*i_vtx  ];
            _entity_center[3*n_elmt+1] += pvtx_coord[3*i_vtx+1];
            _entity_center[3*n_elmt+2] += pvtx_coord[3*i_vtx+2];
          }
          _entity_center[3*n_elmt  ] = _entity_center[3*n_elmt  ] * pond;
          _entity_center[3*n_elmt+1] = _entity_center[3*n_elmt+1] * pond;
          _entity_center[3*n_elmt+2] = _entity_center[3*n_elmt+2] * pond;
          n_elmt++;
        }
      } /* End section */
    }
    shift_part += fctv->n_part[i_domain];
  }

}


static
void
_prepare_vtx_cell_nodal
(
 PDM_field_cell_to_vtx_t*  fctv
)
{

  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;
  int order = 1;
  /*
   * Create vtx_cell
   */
  int shift_part = 0;
  fctv->pvtx_cell_n   = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int *));
  fctv->pvtx_cell_idx = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int *));
  fctv->pvtx_cell     = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int *));
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {

    PDM_part_mesh_nodal_elmts_t* pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(fctv->pmn[i_domain], geom_kind);
    int n_section    = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
    int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {

      int n_vtx        = PDM_part_mesh_nodal_n_vtx_get      (fctv->pmn[i_domain], i_part);

      int n_elmt = 0;
      int n_cell_vtx = 0;
      for(int i_section = 0; i_section < n_section; ++i_section) {
        PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, sections_id[i_section]);
        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
        int n_elmt_in_section = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                                            sections_id[i_section],
                                                                            i_part);
        n_elmt     += n_elmt_in_section;
        n_cell_vtx += n_elmt_in_section * n_vtx_per_elmt;
      }

      assert(n_elmt == fctv->parts[i_domain][i_part].n_cell); // Pour l'instant 3D

      int *_cell_vtx_idx = malloc((n_elmt+1) * sizeof(int));
      int *_cell_vtx     = malloc(n_cell_vtx * sizeof(int));

      n_cell_vtx = 0;
      n_elmt     = 0;
      _cell_vtx_idx[0] = 0;
      for(int i_section = 0; i_section < n_section; ++i_section) {
        PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, sections_id[i_section]);
        int n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
        int n_elmt_in_section = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                                            sections_id[i_section],
                                                                            i_part);
        int            *connec              = NULL;
        PDM_g_num_t    *numabs              = NULL;
        int            *parent_num          = NULL;
        PDM_g_num_t    *parent_entity_g_num = NULL;
        PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                  sections_id[i_section],
                                                  i_part,
                                                  &connec,
                                                  &numabs,
                                                  &parent_num,
                                                  &parent_entity_g_num,
                                                  PDM_OWNERSHIP_KEEP);

        for(int i_elmt = 0; i_elmt < n_elmt_in_section; ++i_elmt) {
          _cell_vtx_idx[n_elmt+1] = _cell_vtx_idx[n_elmt] + n_vtx_per_elmt;
          for(int k = 0; k < n_vtx_per_elmt; ++k) {
            _cell_vtx[n_cell_vtx++] = connec[n_vtx_per_elmt*i_elmt+k];
          }
          n_elmt++;
        }
      } /* End section */

      /* Transpose */
      PDM_connectivity_transpose(n_elmt,
                                 n_vtx,
                                 _cell_vtx_idx,
                                 _cell_vtx,
                                 &fctv->pvtx_cell_idx[i_part+shift_part],
                                 &fctv->pvtx_cell    [i_part+shift_part]);

      fctv->pvtx_cell_n[i_part+shift_part] = malloc(fctv->parts[i_domain][i_part].n_vtx * sizeof(int));
      for(int i_vtx = 0; i_vtx < fctv->parts[i_domain][i_part].n_vtx; ++i_vtx) {
        fctv->pvtx_cell_n[i_part+shift_part][i_vtx] = fctv->pvtx_cell_idx[i_part+shift_part][i_vtx+1] - fctv->pvtx_cell_idx[i_part+shift_part][i_vtx];
      }

      /* Free */
      free(_cell_vtx_idx);
      free(_cell_vtx    );

    }
    shift_part += fctv->n_part[i_domain];
  }

}




static
void
_warm_up_distant_neighbor
(
 PDM_field_cell_to_vtx_t*  fctv
)
{
  /* Deduce graph with all graphe inside same domain and between domain */
  // int ***vtx_part_bound_proc_idx = fctv->entity_part_bound_proc_idx[PDM_MESH_ENTITY_VTX];
  int ***vtx_part_bound_part_idx = fctv->entity_part_bound_part_idx[PDM_MESH_ENTITY_VTX];
  int ***vtx_part_bound          = fctv->entity_part_bound         [PDM_MESH_ENTITY_VTX];


  int         **pdi_neighbor_idx         = NULL;
  int         **pdi_neighbor             = NULL;
  int           n_composed_interface     = 0;
  int          *composed_interface_idx   = NULL;
  int          *composed_interface       = NULL;
  PDM_g_num_t  *composed_ln_to_gn_sorted = NULL;

  if(fctv->pdi != NULL) {

    int is_describe_vtx  = PDM_part_domain_interface_exist_get(fctv->pdi, PDM_BOUND_TYPE_VTX );
    int is_describe_edge = PDM_part_domain_interface_exist_get(fctv->pdi, PDM_BOUND_TYPE_EDGE);
    int is_describe_face = PDM_part_domain_interface_exist_get(fctv->pdi, PDM_BOUND_TYPE_FACE);

    int is_describe_vtx_l  = is_describe_vtx;
    int is_describe_edge_l = is_describe_edge;
    int is_describe_face_l = is_describe_face;
    PDM_MPI_Allreduce(&is_describe_vtx_l , &is_describe_vtx , 1, PDM_MPI_INT, PDM_MPI_MAX, fctv->comm);
    PDM_MPI_Allreduce(&is_describe_edge_l, &is_describe_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, fctv->comm);
    PDM_MPI_Allreduce(&is_describe_face_l, &is_describe_face, 1, PDM_MPI_INT, PDM_MPI_MAX, fctv->comm);

    printf("is_describe_vtx_l  = %i \n", is_describe_vtx_l );
    printf("is_describe_edge_l = %i \n", is_describe_edge_l);
    printf("is_describe_face_l = %i \n", is_describe_face_l);

    int **pn_vtx = malloc(fctv->n_domain * sizeof(int *));
    for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
      pn_vtx[i_domain] = malloc(fctv->n_part[i_domain] * sizeof(int));
      for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
        pn_vtx[i_domain][i_part] = fctv->parts[i_domain][i_part].n_vtx;
      }
    }

    if(is_describe_vtx_l == 0 && is_describe_face_l == 1) {
      int          **pn_face        = malloc(fctv->n_domain * sizeof(int          *));
      PDM_g_num_t ***pface_ln_to_gn = malloc(fctv->n_domain * sizeof(PDM_g_num_t **));
      PDM_g_num_t ***pvtx_ln_to_gn  = malloc(fctv->n_domain * sizeof(PDM_g_num_t **));
      int         ***pface_vtx_idx  = malloc(fctv->n_domain * sizeof(int         **));
      int         ***pface_vtx      = malloc(fctv->n_domain * sizeof(int         **));
      for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
        pn_face       [i_domain] = malloc(fctv->n_part[i_domain] * sizeof(int          ));
        pvtx_ln_to_gn [i_domain] = malloc(fctv->n_part[i_domain] * sizeof(PDM_g_num_t *));
        pface_ln_to_gn[i_domain] = malloc(fctv->n_part[i_domain] * sizeof(PDM_g_num_t *));
        pface_vtx_idx [i_domain] = malloc(fctv->n_part[i_domain] * sizeof(int         *));
        pface_vtx     [i_domain] = malloc(fctv->n_part[i_domain] * sizeof(int         *));
        for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
          pn_face       [i_domain][i_part] = fctv->parts[i_domain][i_part].n_face;
          pvtx_ln_to_gn [i_domain][i_part] = fctv->parts[i_domain][i_part].vtx_ln_to_gn;
          pface_ln_to_gn[i_domain][i_part] = fctv->parts[i_domain][i_part].face_ln_to_gn;
          pface_vtx_idx [i_domain][i_part] = fctv->parts[i_domain][i_part].face_vtx_idx;
          pface_vtx     [i_domain][i_part] = fctv->parts[i_domain][i_part].face_vtx;
        }
      }

      PDM_part_domain_interface_face2vtx(fctv->pdi,
                                         fctv->n_part,
                                         pn_face,
                                         pface_ln_to_gn,
                                         pn_vtx,
                                         pvtx_ln_to_gn,
                                         pface_vtx_idx,
                                         pface_vtx);

      for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
        free(pn_face       [i_domain]);
        free(pvtx_ln_to_gn [i_domain]);
        free(pface_ln_to_gn[i_domain]);
        free(pface_vtx_idx [i_domain]);
        free(pface_vtx     [i_domain]);
      }
      free(pn_face);
      free(pvtx_ln_to_gn );
      free(pface_ln_to_gn);
      free(pface_vtx_idx );
      free(pface_vtx     );
    }


    PDM_part_domain_interface_as_graph(fctv->pdi,
                                       PDM_BOUND_TYPE_VTX,
                                       pn_vtx,
                                       NULL,
                                       &pdi_neighbor_idx,
                                       &pdi_neighbor,
                                       &n_composed_interface,
                                       &composed_interface_idx,
                                       &composed_interface,
                                       &composed_ln_to_gn_sorted);

    for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
      free(pn_vtx        [i_domain]);
    }
    free(pn_vtx);
  }


  int shift_part   = 0;
  int shift_part_g = 0;

  fctv->neighbor_idx       = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int *));
  fctv->neighbor_desc      = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int *));
  fctv->neighbor_interface = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int *));
  fctv->pn_vtx             = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int  ));
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {

    int n_part_total = fctv->n_part_g_idx[i_domain+1] - fctv->n_part_g_idx[i_domain];

    /* First loop to count */
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {

      int n_vtx = fctv->parts[i_domain][i_part].n_vtx;
      fctv->pn_vtx[i_part+shift_part] = n_vtx;

      fctv->neighbor_idx[i_part+shift_part] = malloc((n_vtx+1) * sizeof(int));
      int* _neighbor_idx   = fctv->neighbor_idx[i_part+shift_part];
      int* _vtx_part_bound = vtx_part_bound[i_domain][i_part];

      int* _neighbor_n = PDM_array_zeros_int(n_vtx);

      int n_part_entity_bound_tot = vtx_part_bound_part_idx[i_domain][i_part][n_part_total];
      for(int idx_entity = 0; idx_entity < n_part_entity_bound_tot; ++idx_entity) {
        int i_entity = _vtx_part_bound[4*idx_entity]-1;
        _neighbor_n[i_entity] += 1;
      }

      /*
       * Add comming from interface
       */
      if(pdi_neighbor_idx != NULL) {
        for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
          _neighbor_n[i_entity] += pdi_neighbor_idx[i_part+shift_part][i_entity+1] - pdi_neighbor_idx[i_part+shift_part][i_entity];
        }
      }

      if(0 == 1) {
        PDM_log_trace_array_int(_neighbor_n, n_vtx, "_neighbor_n ::");
      }

      /* Compute index */
      _neighbor_idx[0] = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        _neighbor_idx[i_entity+1] = _neighbor_idx[i_entity] + _neighbor_n[i_entity];
        _neighbor_n[i_entity] = 0;
      }

      if(0 == 1) {
        PDM_log_trace_array_int(_neighbor_idx, n_vtx, "_neighbor_idx ::");
      }

      fctv->neighbor_desc     [i_part+shift_part] = (int *) malloc( 3 * _neighbor_idx[n_vtx] * sizeof(int) );
      fctv->neighbor_interface[i_part+shift_part] = (int *) malloc(     _neighbor_idx[n_vtx] * sizeof(int) );
      int* _neighbor_desc      = fctv->neighbor_desc     [i_part+shift_part];
      int* _neighbor_interface = fctv->neighbor_interface[i_part+shift_part];

      /* Fill */
      for(int idx_entity = 0; idx_entity < vtx_part_bound_part_idx[i_domain][i_part][n_part_total]; ++idx_entity) {
        int i_entity = _vtx_part_bound[4*idx_entity]-1;
        int idx_write = _neighbor_idx[i_entity] + _neighbor_n[i_entity]++;
        _neighbor_desc[3*idx_write  ] = _vtx_part_bound[4*idx_entity+1];                // i_proc_opp;
        _neighbor_desc[3*idx_write+1] = _vtx_part_bound[4*idx_entity+2]+shift_part_g-1; // i_part_opp
        _neighbor_desc[3*idx_write+2] = _vtx_part_bound[4*idx_entity+3]-1;              // i_entity_opp
        _neighbor_interface[idx_write] = -40000;
      }

      if(pdi_neighbor_idx != NULL) {
        for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
          for(int idx = pdi_neighbor_idx[i_part+shift_part][i_entity]; idx < pdi_neighbor_idx[i_part+shift_part][i_entity+1]; ++idx) {
            int idx_write = _neighbor_idx[i_entity] + _neighbor_n[i_entity]++;
            _neighbor_desc[3*idx_write  ]  = pdi_neighbor[i_part+shift_part][4*idx  ];
            _neighbor_desc[3*idx_write+1]  = pdi_neighbor[i_part+shift_part][4*idx+1];
            _neighbor_desc[3*idx_write+2]  = pdi_neighbor[i_part+shift_part][4*idx+2];
            _neighbor_interface[idx_write] = pdi_neighbor[i_part+shift_part][4*idx+3];
          }
        }
      }

      free(_neighbor_n);
    }

    shift_part   += fctv->n_part              [i_domain];
    shift_part_g += n_part_total;
  }

  shift_part = 0;
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {

      if(fctv->pdi != NULL) {
        free(pdi_neighbor_idx[i_part+shift_part]);
        free(pdi_neighbor    [i_part+shift_part]);
      }


    }
    shift_part   += fctv->n_part[i_domain];
  }


  if(fctv->pdi != NULL) {
    free(pdi_neighbor_idx);
    free(pdi_neighbor    );

    free(composed_interface_idx  );
    free(composed_interface      );
    free(composed_ln_to_gn_sorted);
  }

}

static
void
_create_bnd_graph
(
 PDM_field_cell_to_vtx_t* fctv
)
{
  /*
   * Create also vtx_face_bound :
   *   - Vertex can be on 2 boundaries
   */
  int shift_part = 0;
  fctv->vtx_face_bound_idx    = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));
  fctv->vtx_face_bound_n      = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));
  fctv->vtx_face_bound        = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));
  fctv->vtx_face_bound_group  = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));
  fctv->vtx_face_bound_coords = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));

  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {

    int n_group = fctv->n_group[i_domain];

    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {

      int n_vtx = fctv->parts[i_domain][i_part].n_vtx;
      fctv->vtx_face_bound_idx[i_part+shift_part] = malloc( (n_vtx + 1) * sizeof(int));
      fctv->vtx_face_bound_n  [i_part+shift_part] = PDM_array_zeros_int(n_vtx);

      int *_vtx_face_bound_idx = (int *) fctv->vtx_face_bound_idx[i_part+shift_part];
      int *_vtx_face_bound_n   = (int *) fctv->vtx_face_bound_n  [i_part+shift_part];

      int  *n_face_group = fctv->n_group_entity[PDM_BOUND_TYPE_FACE][i_domain][i_part];
      int **face_group   = fctv->group_entity  [PDM_BOUND_TYPE_FACE][i_domain][i_part];

      int* face_vtx_idx = fctv->parts[i_domain][i_part].face_vtx_idx;
      int* face_vtx     = fctv->parts[i_domain][i_part].face_vtx;

      for(int i_group = 0; i_group < n_group; ++i_group) {
        for(int idx_face = 0; idx_face < n_face_group[i_group]; ++idx_face) {
          int i_face = face_group[i_group][idx_face]-1;
          for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = face_vtx[idx_vtx]-1;
            _vtx_face_bound_n[i_vtx] += 1;
          }
        }
      }

      _vtx_face_bound_idx[0] = 0;
      for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
        _vtx_face_bound_idx[i_vtx+1] = _vtx_face_bound_idx[i_vtx] + _vtx_face_bound_n[i_vtx];
        _vtx_face_bound_n[i_vtx] = 0;
      }
      fctv->vtx_face_bound       [i_part+shift_part] = malloc(    _vtx_face_bound_idx[n_vtx] * sizeof(int   ));
      fctv->vtx_face_bound_group [i_part+shift_part] = malloc(    _vtx_face_bound_idx[n_vtx] * sizeof(int   ));
      fctv->vtx_face_bound_coords[i_part+shift_part] = malloc(3 * _vtx_face_bound_idx[n_vtx] * sizeof(double));
      int    *_vtx_face_bound        = fctv->vtx_face_bound       [i_part+shift_part];
      int    *_vtx_face_bound_group  = fctv->vtx_face_bound_group [i_part+shift_part];
      double *_vtx_face_bound_coords = fctv->vtx_face_bound_coords[i_part+shift_part];
      double *_pvtx_coord            = fctv->parts[i_domain][i_part].vtx;

      for(int i_group = 0; i_group < n_group; ++i_group) {
        for(int idx_face = 0; idx_face < n_face_group[i_group]; ++idx_face) {
          int i_face = face_group[i_group][idx_face]-1;

          double center_face[3] = {0., 0., 0.};
          double pond = 1./((double)  face_vtx_idx[i_face+1] - face_vtx_idx[i_face]);
          for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = face_vtx[idx_vtx]-1;
            for(int k = 0; k < 3; ++k) {
              center_face[k] += _pvtx_coord[3*i_vtx+k];
            }
          }
          for(int k = 0; k < 3; ++k) {
            center_face[k] = center_face[k] * pond;
          }

          for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
            int i_vtx = face_vtx[idx_vtx]-1;
            int idx_write = _vtx_face_bound_idx[i_vtx] + _vtx_face_bound_n[i_vtx]++;
            _vtx_face_bound      [idx_write] = idx_face;
            _vtx_face_bound_group[idx_write] = i_group;
            for(int k = 0; k < 3; ++k) {
              _vtx_face_bound_coords[3*idx_write+k] = center_face[k];
            }
          }
        }
      }
    }
    shift_part   += fctv->n_part              [i_domain];
  }
}

static
void
_create_bnd_graph_nodal
(
 PDM_field_cell_to_vtx_t* fctv
)
{

  /*
   * Create also vtx_face_bound :
   *   - Vertex can be on 2 boundaries
   */
  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_SURFACIC;
  int order = 1;

  int shift_part = 0;
  fctv->vtx_face_bound_idx    = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));
  fctv->vtx_face_bound_n      = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));
  fctv->vtx_face_bound        = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));
  fctv->vtx_face_bound_group  = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));
  fctv->vtx_face_bound_coords = malloc(fctv->n_part_g_idx[fctv->n_domain] * sizeof(int * ));

  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {

    PDM_part_mesh_nodal_elmts_t* pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(fctv->pmn[i_domain], geom_kind);

    int n_section    = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
    int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

    int **connect        = malloc(n_section * sizeof(int *));
    int  *n_vtx_per_elmt = malloc(n_section * sizeof(int  ));

    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {

      int n_vtx = fctv->parts[i_domain][i_part].n_vtx;

      /* Setup idx */
      int *section_elmt_idx = PDM_part_mesh_nodal_elmts_compute_sections_idx(pmne,
                                                                       i_part);

      int n_group_part = PDM_part_mesh_nodal_elmts_n_group_get(pmne);


      fctv->vtx_face_bound_idx[i_part+shift_part] = malloc( (n_vtx + 1) * sizeof(int));
      fctv->vtx_face_bound_n  [i_part+shift_part] = PDM_array_zeros_int(n_vtx);

      int *_vtx_face_bound_idx = (int *) fctv->vtx_face_bound_idx[i_part+shift_part];
      int *_vtx_face_bound_n   = (int *) fctv->vtx_face_bound_n  [i_part+shift_part];

      /*
       * Prepare usefull shorcut
       */
      for(int i_section = 0; i_section < n_section; ++i_section) {
        PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, sections_id[i_section]);
        n_vtx_per_elmt[i_section] = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

        PDM_g_num_t    *numabs              = NULL;
        int            *parent_num          = NULL;
        PDM_g_num_t    *parent_entity_g_num = NULL;
        PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                  sections_id[i_section],
                                                  i_part,
                                                  &connect[i_section],
                                                  &numabs,
                                                  &parent_num,
                                                  &parent_entity_g_num,
                                                  PDM_OWNERSHIP_KEEP);
      }

      for(int i_group = 0; i_group < n_group_part; ++i_group) {

        int          n_group_elmt   = 0;
        int         *group_elmt     = NULL;
        PDM_g_num_t *group_ln_to_gn = NULL;

        PDM_part_mesh_nodal_elmts_group_get(pmne,
                                            i_part,
                                            i_group,
                                            &n_group_elmt,
                                            &group_elmt,
                                            &group_ln_to_gn,
                                            PDM_OWNERSHIP_KEEP);

        for(int idx_elmt = 0; idx_elmt < n_group_elmt; ++idx_elmt) {
          int i_elmt    = group_elmt[idx_elmt]-1;
          int i_section = PDM_binary_search_gap_int(i_elmt, section_elmt_idx, n_section+1);
          int i_loc_elmt = i_elmt - section_elmt_idx[i_section];
          int ln_vtx_per_elmt = n_vtx_per_elmt[i_section];
          for(int k = 0; k < ln_vtx_per_elmt; ++k) {
            int i_vtx = connect[i_section][ln_vtx_per_elmt*i_loc_elmt+k]-1;
            _vtx_face_bound_n[i_vtx] += 1;
          }
        }
      }

      _vtx_face_bound_idx[0] = 0;
      for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
        _vtx_face_bound_idx[i_vtx+1] = _vtx_face_bound_idx[i_vtx] + _vtx_face_bound_n[i_vtx];
        _vtx_face_bound_n[i_vtx] = 0;
      }

      fctv->vtx_face_bound       [i_part+shift_part] = malloc(    _vtx_face_bound_idx[n_vtx] * sizeof(int   ));
      fctv->vtx_face_bound_group [i_part+shift_part] = malloc(    _vtx_face_bound_idx[n_vtx] * sizeof(int   ));
      fctv->vtx_face_bound_coords[i_part+shift_part] = malloc(3 * _vtx_face_bound_idx[n_vtx] * sizeof(double));
      int    *_vtx_face_bound        = fctv->vtx_face_bound       [i_part+shift_part];
      int    *_vtx_face_bound_group  = fctv->vtx_face_bound_group [i_part+shift_part];
      double *_vtx_face_bound_coords = fctv->vtx_face_bound_coords[i_part+shift_part];
      double *_pvtx_coord            = fctv->parts[i_domain][i_part].vtx;

      for(int i_group = 0; i_group < n_group_part; ++i_group) {

        int          n_group_elmt   = 0;
        int         *group_elmt     = NULL;
        PDM_g_num_t *group_ln_to_gn = NULL;

        PDM_part_mesh_nodal_elmts_group_get(pmne,
                                            i_part,
                                            i_group,
                                            &n_group_elmt,
                                            &group_elmt,
                                            &group_ln_to_gn,
                                            PDM_OWNERSHIP_KEEP);

        for(int idx_elmt = 0; idx_elmt < n_group_elmt; ++idx_elmt) {
          int i_elmt    = group_elmt[idx_elmt]-1;
          int i_section = PDM_binary_search_gap_int(i_elmt, section_elmt_idx, n_section+1);
          int i_loc_elmt = i_elmt - section_elmt_idx[i_section];
          int ln_vtx_per_elmt = n_vtx_per_elmt[i_section];

          double center_face[3] = {0., 0., 0.};
          double pond = 1./((double) ln_vtx_per_elmt);

          for(int k = 0; k < ln_vtx_per_elmt; ++k) {
            int i_vtx = connect[i_section][ln_vtx_per_elmt*i_loc_elmt+k]-1;
            for(int p = 0; p < 3; ++p) {
              center_face[p] += _pvtx_coord[3*i_vtx+p];
            }
          }
          for(int k = 0; k < 3; ++k) {
            center_face[k] = center_face[k] * pond;
          }

          for(int k = 0; k < ln_vtx_per_elmt; ++k) {
            int i_vtx = connect[i_section][ln_vtx_per_elmt*i_loc_elmt+k]-1;
            int idx_write = _vtx_face_bound_idx[i_vtx] + _vtx_face_bound_n[i_vtx]++;
            _vtx_face_bound      [idx_write] = idx_elmt; // - group_elmt_idx[i_group];
            _vtx_face_bound_group[idx_write] = i_group;
            for(int p = 0; p < 3; ++p) {
              _vtx_face_bound_coords[3*idx_write+p] = center_face[p];
            }
          }
        }
      }

      free(section_elmt_idx);
    }

    free(connect       );
    free(n_vtx_per_elmt);
    shift_part   += fctv->n_part              [i_domain];
  }

}


static
inline
double
_evaluate_distance
(
  double x1[3],
  double x2[3],
  int    p
)
{
  double dp = (double)  p;
  double dist = 0.;
  for(int k = 0; k < 3; ++k) {
    dist += pow(x2[k] - x1[k], dp);
  }
  if(p > 0){
    double ip = 1./p;
    dist = pow(dist, ip);
  }
  return  dist;
}


static
void
_interpolate_one_part
(
  int      n_vtx,
  double  *vtx_coord,
  double  *cell_center,
  int      stride,
  int     *vtx_cell_idx,
  int     *vtx_cell,
  double  *plocal_field,
  double **pbound_field,
  int     *neighbor_idx,
  int     *vtx_face_bound_idx,
  int     *vtx_face_bound,
  int     *vtx_face_bound_group,
  double  *vtx_face_bound_coords,
  int     *pvtx_cell_coords_opp_n,
  double  *pvtx_cell_coords_opp,
  int     *pvtx_face_bound_coords_opp_n,
  double  *pvtx_face_bound_coords_opp,
  int     *pvtx_cell_field_opp_n,
  double  *pvtx_cell_field_opp,
  int     *pvtx_face_field_opp_n,
  double  *pvtx_face_field_opp,
  double **result_field
)
{
  PDM_UNUSED(n_vtx);
  PDM_UNUSED(stride);
  PDM_UNUSED(vtx_cell_idx);
  PDM_UNUSED(vtx_cell);
  PDM_UNUSED(plocal_field);
  PDM_UNUSED(pbound_field);
  PDM_UNUSED(neighbor_idx);
  PDM_UNUSED(pvtx_cell_coords_opp_n);
  PDM_UNUSED(pvtx_cell_coords_opp);
  PDM_UNUSED(pvtx_face_bound_coords_opp_n);
  PDM_UNUSED(pvtx_face_bound_coords_opp);
  PDM_UNUSED(pvtx_cell_field_opp_n);
  PDM_UNUSED(pvtx_cell_field_opp);
  PDM_UNUSED(pvtx_face_field_opp_n);
  PDM_UNUSED(pvtx_face_field_opp);

  int *pvtx_cell_coords_opp_idx       = PDM_array_new_idx_from_sizes_int(pvtx_cell_coords_opp_n      , neighbor_idx[n_vtx]);
  int *pvtx_face_bound_coords_opp_idx = PDM_array_new_idx_from_sizes_int(pvtx_face_bound_coords_opp_n, neighbor_idx[n_vtx]);

  int p = 0;

  double *_result_field = malloc(stride * n_vtx * sizeof(double));

  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx) {

    // double vtx_coord_x = vtx_coord[3*i_vtx  ];
    // double vtx_coord_y = vtx_coord[3*i_vtx+1];
    // double vtx_coord_z = vtx_coord[3*i_vtx+2];

    for(int k = 0; k < stride; ++k) {
      _result_field[stride*i_vtx+k] =  0.;
    }
    double tot_dist = 0.;

    /* Pour l'instant ->  Inverse distance weighting */

    /* Interior */
    for(int idx_cell = vtx_cell_idx[i_vtx]; idx_cell < vtx_cell_idx[i_vtx+1]; ++idx_cell) {
      int i_cell = PDM_ABS(vtx_cell[idx_cell])-1;

      // double dx = cell_center[3*i_cell  ] - vtx_coord_x;
      // double dy = cell_center[3*i_cell+1] - vtx_coord_y;
      // double dz = cell_center[3*i_cell+2] - vtx_coord_z;
      // double dist = sqrt(dx * dx + dy * dy + dz * dz);
      double dist = _evaluate_distance(&cell_center[3*i_cell], &vtx_coord[3*i_vtx], p);

      // printf("dist  = %12.5e / plocal = %12.5e\n", dist, plocal_field[i_cell]);

      double inv_dist = 1./PDM_MAX(dist, 1.e-12);
      tot_dist += inv_dist;

      for(int k = 0; k < stride; ++k) {
        _result_field[stride*i_vtx+k] += inv_dist * plocal_field[stride*i_cell+k];
      }
    }

    /* Interior BND */
    for(int idx_face = vtx_face_bound_idx[i_vtx]; idx_face < vtx_face_bound_idx[i_vtx+1]; ++idx_face) {
      int i_face  = PDM_ABS(vtx_face_bound      [idx_face]);
      int i_group = PDM_ABS(vtx_face_bound_group[idx_face]);

      // double dx = vtx_face_bound_coords[3*idx_face  ] - vtx_coord_x;
      // double dy = vtx_face_bound_coords[3*idx_face+1] - vtx_coord_y;
      // double dz = vtx_face_bound_coords[3*idx_face+2] - vtx_coord_z;
      // double dist = sqrt(dx * dx + dy * dy + dz * dz);
      double dist = _evaluate_distance(&vtx_face_bound_coords[3*idx_face  ], &vtx_coord[3*i_vtx], p);

      double inv_dist = 1./PDM_MAX(dist, 1.e-12);
      tot_dist += inv_dist;

      for(int k = 0; k < stride; ++k) {
        _result_field[stride*i_vtx+k] += inv_dist * pbound_field[i_group][stride*i_face+k];
      }
    }

    /* Distant cell */
    for(int idx_neigh = neighbor_idx[i_vtx]; idx_neigh < neighbor_idx[i_vtx+1];  ++idx_neigh)  {

      assert(pvtx_cell_field_opp_n[idx_neigh] == pvtx_cell_field_opp_n[idx_neigh]);
       for(int idx_recv = pvtx_cell_coords_opp_idx[idx_neigh]; idx_recv < pvtx_cell_coords_opp_idx[idx_neigh+1]; ++idx_recv){

        // double dx = pvtx_cell_coords_opp[3*idx_recv  ] - vtx_coord_x;
        // double dy = pvtx_cell_coords_opp[3*idx_recv+1] - vtx_coord_y;
        // double dz = pvtx_cell_coords_opp[3*idx_recv+2] - vtx_coord_z;
        // double dist = sqrt(dx * dx + dy * dy + dz * dz);
        double dist = _evaluate_distance(&pvtx_cell_coords_opp[3*idx_recv  ], &vtx_coord[3*i_vtx], p);

        double inv_dist = 1./PDM_MAX(dist, 1.e-12);
        tot_dist += inv_dist;

        for(int k = 0; k < stride; ++k) {
          _result_field[stride*i_vtx+k] += inv_dist * pvtx_cell_field_opp[stride*idx_recv+k];
        }
      }
    }

    /* Distant BND */
    for(int idx_neigh = neighbor_idx[i_vtx]; idx_neigh < neighbor_idx[i_vtx+1];  ++idx_neigh)  {

      assert(pvtx_face_bound_coords_opp_n[idx_neigh] == pvtx_face_field_opp_n[idx_neigh]);
       for(int idx_recv = pvtx_face_bound_coords_opp_idx[idx_neigh]; idx_recv < pvtx_face_bound_coords_opp_idx[idx_neigh+1]; ++idx_recv){

        // double dx = pvtx_face_bound_coords_opp[3*idx_recv  ] - vtx_coord_x;
        // double dy = pvtx_face_bound_coords_opp[3*idx_recv+1] - vtx_coord_y;
        // double dz = pvtx_face_bound_coords_opp[3*idx_recv+2] - vtx_coord_z;
        // double dist = sqrt(dx * dx + dy * dy + dz * dz);
        double dist = _evaluate_distance(&pvtx_face_bound_coords_opp[3*idx_recv  ], &vtx_coord[3*i_vtx], p);

        double inv_dist = 1./PDM_MAX(dist, 1.e-12);
        tot_dist += inv_dist;

        for(int k = 0; k < stride; ++k) {
          _result_field[stride*i_vtx+k] += inv_dist * pvtx_face_field_opp[stride*idx_recv+k];
        }
      }
    }

    /* Finalize ponderate */
    double inv_tot_dist = 1./tot_dist;
    for(int k = 0; k < stride; ++k) {
      _result_field[stride*i_vtx+k] = _result_field[stride*i_vtx+k] * inv_tot_dist;
    }

  }

  *result_field = _result_field;
  free(pvtx_cell_coords_opp_idx      );
  free(pvtx_face_bound_coords_opp_idx);
}

static
int
_field_kind_to_size
(
 PDM_field_kind_t  field_kind
)
{
  if(field_kind  == PDM_FIELD_KIND_SCALAR) {
    return 1;
  } else  if(field_kind  == PDM_FIELD_KIND_COORDS) {
    return 3;
  } else  if(field_kind  == PDM_FIELD_KIND_VECTOR) {
    return 3;
  } else  if(field_kind  == PDM_FIELD_KIND_TENSOR_SYM) {
    return 6;
  }
  return -1;
}


static
void
_apply_transformation
(
  PDM_part_domain_interface_t  *pdi,
  PDM_field_kind_t              field_kind,
  int                           n_vtx,
  int                          *neighbor_idx,
  int                          *neighbor_interface,
  int                          *pvtx_cell_field_opp_n,
  double                       *pvtx_cell_field_opp
)
{
  int nrecv = 0;
  double rot[3][3];

  if(field_kind == PDM_FIELD_KIND_COORDS) {
    // assert(stride ==  1);
    /* Apply translation AND  Rotation if any */
    for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
      for(int idx_entity = neighbor_idx[i_entity]; idx_entity < neighbor_idx[i_entity+1]; ++idx_entity) {

        if(neighbor_interface[idx_entity] != -40000) {
          int  i_interface = PDM_ABS(neighbor_interface[idx_entity])-1;

          if(pdi->translation_vect[i_interface] != NULL) {
            for(int idx_recv = 0; idx_recv < pvtx_cell_field_opp_n[idx_entity]; ++idx_recv){
              for(int k = 0; k < 3; ++k) {
                pvtx_cell_field_opp[3*(nrecv+idx_recv)+k] += PDM_SIGN(neighbor_interface[idx_entity]) * pdi->translation_vect[i_interface][k];
              }
            }
          } else if(pdi->rotation_direction[i_interface] != NULL){

            for(int idx_recv = 0; idx_recv < pvtx_cell_field_opp_n[idx_entity]; ++idx_recv){
              double x = pvtx_cell_field_opp[3*(nrecv+idx_recv)  ];
              double y = pvtx_cell_field_opp[3*(nrecv+idx_recv)+1];
              double z = pvtx_cell_field_opp[3*(nrecv+idx_recv)+2];

              double xp = x - pdi->rotation_center[i_interface][0];
              double yp = y - pdi->rotation_center[i_interface][1];
              double zp = z - pdi->rotation_center[i_interface][2];

              double c   = cos(PDM_SIGN(neighbor_interface[idx_entity]) * pdi->rotation_angle[i_interface]);
              double s   = sin(PDM_SIGN(neighbor_interface[idx_entity]) * pdi->rotation_angle[i_interface]);
              double omc = 1.-c;
              double xa  = pdi->rotation_direction[i_interface][0];
              double ya  = pdi->rotation_direction[i_interface][1];
              double za  = pdi->rotation_direction[i_interface][2];

              rot[0][0] = xa * xa * omc + c;
              rot[0][1] = xa * ya * omc - za * s;
              rot[0][2] = xa * za * omc + ya * s;

              rot[1][0] = xa * ya * omc + za * s;
              rot[1][1] = ya * ya * omc + c;
              rot[1][2] = ya * za * omc - xa * s;

              rot[2][0] = xa * za * omc - ya * s;
              rot[2][1] = ya * za * omc + xa * s;
              rot[2][2] = za * za * omc + c;

              for(int k = 0; k < 3; ++k) {
                pvtx_cell_field_opp[3*(nrecv+idx_recv)+k] = pdi->rotation_center[i_interface][k] + (rot[k][0]*xp + rot[k][1]*yp + rot[k][2]*zp);
              }
            } /*  End idx_recv */
          }
        }
        nrecv += pvtx_cell_field_opp_n[idx_entity];
      }
    }
  } else if(field_kind == PDM_FIELD_KIND_VECTOR)  {





  }


}

static
void
_interpolate
(
  PDM_field_cell_to_vtx_t    *fctv,
  PDM_field_kind_t            field_kind,
  int                         stride,
  double                   ***plocal_field,
  double                  ****pbound_field,
  int                       **pvtx_cell_field_opp_n,
  double                    **pvtx_cell_field_opp,
  int                       **pvtx_face_field_opp_n,
  double                    **pvtx_face_field_opp,
  double                  ****result_field
)
{

  double ***_result_field = malloc(fctv->n_domain * sizeof(double **));
  int shift_part = 0;
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    /* First loop to count */
    _result_field[i_domain] = malloc(fctv->n_part[i_domain] * sizeof(double *));
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {

      /* Apply transformation  */
      _apply_transformation(fctv->pdi,
                            field_kind,
                            fctv->pn_vtx                      [i_part+shift_part],
                            fctv->neighbor_idx                [i_part+shift_part],
                            fctv->neighbor_interface          [i_part+shift_part],
                            pvtx_cell_field_opp_n           [i_part+shift_part],
                            pvtx_cell_field_opp             [i_part+shift_part]);

      _apply_transformation(fctv->pdi,
                            field_kind,
                            fctv->pn_vtx                      [i_part+shift_part],
                            fctv->neighbor_idx                [i_part+shift_part],
                            fctv->neighbor_interface          [i_part+shift_part],
                            pvtx_face_field_opp_n           [i_part+shift_part],
                            pvtx_face_field_opp             [i_part+shift_part]);

      _interpolate_one_part(fctv->pn_vtx                      [i_part+shift_part],
                            fctv->parts                       [i_domain][i_part].vtx,
                            fctv->cell_center                 [i_part+shift_part],
                            stride,
                            fctv->pvtx_cell_idx               [i_part+shift_part],
                            fctv->pvtx_cell                   [i_part+shift_part],
                            plocal_field                    [i_domain][i_part],
                            pbound_field                    [i_domain][i_part],
                            fctv->neighbor_idx                [i_part+shift_part],
                            fctv->vtx_face_bound_idx          [i_part+shift_part],
                            fctv->vtx_face_bound              [i_part+shift_part],
                            fctv->vtx_face_bound_group        [i_part+shift_part],
                            fctv->vtx_face_bound_coords       [i_part+shift_part],
                            fctv->pvtx_cell_coords_opp_n      [i_part+shift_part],
                            fctv->pvtx_cell_coords_opp        [i_part+shift_part],
                            fctv->pvtx_face_bound_coords_opp_n[i_part+shift_part],
                            fctv->pvtx_face_bound_coords_opp  [i_part+shift_part],
                            pvtx_cell_field_opp_n           [i_part+shift_part],
                            pvtx_cell_field_opp             [i_part+shift_part],
                            pvtx_face_field_opp_n           [i_part+shift_part],
                            pvtx_face_field_opp             [i_part+shift_part],
                            &_result_field[i_domain][i_part]);

    }
    shift_part += fctv->n_part[i_domain];
  }


  *result_field = _result_field;

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that compute a global mean
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to \ref PDM_field_cell_to_vtx object
 */
PDM_field_cell_to_vtx_t*
PDM_field_cell_to_vtx_create
(
 const int                            n_domain,
 const int                           *n_part,
 const int                           *n_group,
 const PDM_cell_to_vtx_interp_kind_t  interp_kind, // IDW(p), RBF, LSQ, USER
 const int                            n_depth,
 const PDM_MPI_Comm                   comm
)
{
  PDM_field_cell_to_vtx_t* fctv = malloc(sizeof(PDM_field_cell_to_vtx_t));

  fctv->n_depth     = n_depth;
  fctv->interp_kind = interp_kind;
  fctv->idw_p       = 0;

  fctv->is_nodal    = 0;

  fctv->n_domain    = n_domain;
  fctv->n_part      = malloc( n_domain * sizeof(int)); // Make a copy to avoid pb in cython
  fctv->n_group     = malloc( n_domain * sizeof(int)); // Make a copy to avoid pb in cython
  for(int i = 0; i < fctv->n_domain; ++i) {
    fctv->n_part [i] = n_part [i];
    fctv->n_group[i] = n_group[i];
  }
  fctv->comm        = comm;

  fctv->n_part_idx    = (int * ) malloc( (n_domain + 1) * sizeof(int));
  fctv->n_part_g_idx  = (int * ) malloc( (n_domain + 1) * sizeof(int));
  fctv->parts = malloc(n_domain * sizeof(_part_t *));
  fctv->pmn   = malloc(n_domain * sizeof(PDM_part_mesh_nodal_t *));
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    fctv->pmn[i_domain] = NULL;
  }


  fctv->n_part_idx  [0] = 0;
  fctv->n_part_g_idx[0] = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    fctv->parts[i_domain] = malloc( n_part[i_domain] * sizeof(_part_t));
    fctv->n_part_idx[i_domain+1] = fctv->n_part_idx[i_domain] + fctv->n_part[i_domain];

    int n_part_l = n_part[i_domain];
    int n_part_g = -100;
    PDM_MPI_Allreduce(&n_part_l, &n_part_g, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    fctv->n_part_g_idx[i_domain+1] = fctv->n_part_idx[i_domain] + n_part_g;

  }

  /* Si multidomain on fait un shift et tt roule */
  fctv->n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    fctv->n_part_loc_all_domain += fctv->n_part[i_domain];
  }

  fctv->pn_vtx              = NULL;
  fctv->pvtx_cell_n         = NULL;
  fctv->pvtx_cell_idx       = NULL;
  fctv->pvtx_cell           = NULL;

  fctv->neighbor_idx       = NULL;
  fctv->neighbor_desc      = NULL;
  fctv->neighbor_interface = NULL;

  fctv->vtx_face_bound_idx    = NULL;
  fctv->vtx_face_bound_n      = NULL;
  fctv->vtx_face_bound        = NULL;
  fctv->vtx_face_bound_group  = NULL;
  fctv->vtx_face_bound_coords = NULL;

  fctv->pvtx_cell_coords_opp_n       = NULL;
  fctv->pvtx_cell_coords_opp         = NULL;
  fctv->pvtx_face_bound_coords_opp_n = NULL;
  fctv->pvtx_face_bound_coords_opp   = NULL;

  fctv->dn = NULL;

  /* Graph comm  */
  fctv->entity_part_bound_proc_idx = malloc( PDM_MESH_ENTITY_MAX * sizeof(int ***));
  fctv->entity_part_bound_part_idx = malloc( PDM_MESH_ENTITY_MAX * sizeof(int ***));
  fctv->entity_part_bound          = malloc( PDM_MESH_ENTITY_MAX * sizeof(int ***));
  fctv->graph_comm_is_defined      = malloc( PDM_MESH_ENTITY_MAX * sizeof(int *  ));

  for(int i_kind = 0; i_kind < PDM_MESH_ENTITY_MAX; ++i_kind) {
    fctv->entity_part_bound_proc_idx[i_kind] = malloc(n_domain * sizeof(int **));
    fctv->entity_part_bound_part_idx[i_kind] = malloc(n_domain * sizeof(int **));
    fctv->entity_part_bound         [i_kind] = malloc(n_domain * sizeof(int **));
    fctv->graph_comm_is_defined     [i_kind] = 0;
    for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
      fctv->entity_part_bound_proc_idx[i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      fctv->entity_part_bound_part_idx[i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      fctv->entity_part_bound         [i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
        fctv->entity_part_bound_proc_idx[i_kind][i_domain][i_part] = NULL;
        fctv->entity_part_bound_part_idx[i_kind][i_domain][i_part] = NULL;
        fctv->entity_part_bound         [i_kind][i_domain][i_part] = NULL;
      }
    }
  }

  fctv->pdi = NULL;


  fctv->n_group_entity   = malloc( PDM_GEOMETRY_KIND_MAX * sizeof(int  ***));
  fctv->group_entity     = malloc( PDM_GEOMETRY_KIND_MAX * sizeof(int ****));
  fctv->group_is_defined = malloc( PDM_GEOMETRY_KIND_MAX * sizeof(int *   ));

  for(int i_kind = 0; i_kind < PDM_GEOMETRY_KIND_MAX; ++i_kind) {
    fctv->n_group_entity  [i_kind] = malloc(n_domain * sizeof(int  **));
    fctv->group_entity    [i_kind] = malloc(n_domain * sizeof(int ***));
    fctv->group_is_defined[i_kind] = 0;
    for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
      fctv->n_group_entity  [i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int  *));
      fctv->group_entity    [i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int **));
      for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
        fctv->n_group_entity  [i_kind][i_domain][i_part] = malloc(n_group[i_domain] * sizeof(int  ));
        fctv->group_entity    [i_kind][i_domain][i_part] = malloc(n_group[i_domain] * sizeof(int *));
      }
    }
  }

  fctv->cell_center =  NULL;

  return fctv;
}

void
PDM_field_cell_to_vtx_compute
(
  PDM_field_cell_to_vtx_t *fctv
)
{

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(fctv->comm, &i_rank);
  PDM_MPI_Comm_size(fctv->comm, &n_rank);

  /*
   * Compute graph comm from gnum with PDM_part_generate_entity_graph_comm is not provided
   */
  int         ***pvtx_proc_bound_idx  = NULL;
  int         ***pvtx_part_bound_idx  = NULL;
  int         ***pvtx_bound           = NULL;
  int         ***pvtx_priority        = NULL;
  if(fctv->graph_comm_is_defined[PDM_MESH_ENTITY_VTX] == 0) {
    pvtx_proc_bound_idx = malloc( fctv->n_domain * sizeof(int **));
    pvtx_part_bound_idx = malloc( fctv->n_domain * sizeof(int **));
    pvtx_bound          = malloc( fctv->n_domain * sizeof(int **));
    pvtx_priority       = malloc( fctv->n_domain * sizeof(int **));
    for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {

      int          *pn_vtx        = malloc(fctv->n_part[i_domain] * sizeof(int          ));
      PDM_g_num_t **pvtx_ln_to_gn = malloc(fctv->n_part[i_domain] * sizeof(PDM_g_num_t *));

      for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
        pn_vtx       [i_part] = fctv->parts[i_domain][i_part].n_vtx;
        pvtx_ln_to_gn[i_part] = fctv->parts[i_domain][i_part].vtx_ln_to_gn;
      }

      PDM_g_num_t *part_distribution = PDM_compute_entity_distribution(fctv->comm, fctv->n_part[i_domain] );
      PDM_part_generate_entity_graph_comm(fctv->comm,
                                          part_distribution,
                                          NULL,
                                          fctv->n_part[i_domain],
                                          pn_vtx,
                   (const PDM_g_num_t **) pvtx_ln_to_gn,
                                          NULL,
                                          &pvtx_proc_bound_idx[i_domain],
                                          &pvtx_part_bound_idx[i_domain],
                                          &pvtx_bound         [i_domain],
                                          &pvtx_priority      [i_domain]);
      free(part_distribution);
      free(pn_vtx       );
      free(pvtx_ln_to_gn);

      for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
        fctv->entity_part_bound_proc_idx[PDM_MESH_ENTITY_VTX][i_domain][i_part] = pvtx_proc_bound_idx[i_domain][i_part];
        fctv->entity_part_bound_part_idx[PDM_MESH_ENTITY_VTX][i_domain][i_part] = pvtx_part_bound_idx[i_domain][i_part];
        fctv->entity_part_bound         [PDM_MESH_ENTITY_VTX][i_domain][i_part] = pvtx_bound         [i_domain][i_part];
      }
    }
  }

  if(fctv->is_nodal == 0) {
    _prepare_cell_center     (fctv);
    _prepare_vtx_cell        (fctv);
    _create_bnd_graph        (fctv);
  } else {
    _prepare_cell_center_nodal(fctv);
    _prepare_vtx_cell_nodal   (fctv);
    _create_bnd_graph_nodal   (fctv);
  }
  _warm_up_distant_neighbor(fctv);

  /* Free unused */
  if(fctv->graph_comm_is_defined[PDM_MESH_ENTITY_VTX] == 0) {

    for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
        free(pvtx_proc_bound_idx[i_domain][i_part]);
        free(pvtx_part_bound_idx[i_domain][i_part]);
        free(pvtx_bound         [i_domain][i_part]);
        free(pvtx_priority      [i_domain][i_part]);
      }
      free(pvtx_proc_bound_idx[i_domain]);
      free(pvtx_part_bound_idx[i_domain]);
      free(pvtx_bound         [i_domain]);
      free(pvtx_priority      [i_domain]);
    }
    free(pvtx_proc_bound_idx);
    free(pvtx_part_bound_idx);
    free(pvtx_bound         );
    free(pvtx_priority      );
  }

  int **pvtx_cell_n   = fctv->pvtx_cell_n;
  int **pvtx_cell_idx = fctv->pvtx_cell_idx;
  int **pvtx_cell     = fctv->pvtx_cell;

  /*
   * Create distant_neighbor
   */
  assert(fctv->dn == NULL);
  fctv->dn = PDM_distant_neighbor_create(fctv->comm,
                                         fctv->n_part_loc_all_domain,
                                         fctv->pn_vtx,
                                         fctv->neighbor_idx,
                                         fctv->neighbor_desc);

  /* Prepare coordinates to send */
  int    **pvtx_cell_coords_n =  malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  double **pvtx_cell_coords   =  malloc(fctv->n_part_loc_all_domain * sizeof(double *));
  int shift_part = 0;
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    /* First loop to count */
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
      int n_vtx          = fctv->pn_vtx      [i_part+shift_part];
      int* _neighbor_idx = fctv->neighbor_idx[i_part+shift_part];

      pvtx_cell_coords_n[i_part+shift_part] = PDM_array_zeros_int(n_vtx);
      int *_pvtx_cell_coords_n = pvtx_cell_coords_n[i_part+shift_part];
      int *_pvtx_cell_n        = pvtx_cell_n       [i_part+shift_part];
      int *_pvtx_cell_idx      = pvtx_cell_idx     [i_part+shift_part];
      int *_pvtx_cell          = pvtx_cell         [i_part+shift_part];

      double*  _pcell_center = fctv->cell_center[i_part+shift_part];

      /* Count */
      int n_vtx_cell_to_send = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        int n_neight = _neighbor_idx[i_entity+1] - _neighbor_idx[i_entity];
        if(n_neight > 0 ) {
          n_vtx_cell_to_send += _pvtx_cell_n[i_entity];
          _pvtx_cell_coords_n[i_entity] = _pvtx_cell_n[i_entity];
        }
      }

      pvtx_cell_coords[i_part+shift_part] = malloc( 3 * n_vtx_cell_to_send * sizeof(double));
      double *_pvtx_cell_coords = pvtx_cell_coords[i_part+shift_part];

      n_vtx_cell_to_send = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        int n_neight = _neighbor_idx[i_entity+1] - _neighbor_idx[i_entity];
        if(n_neight > 0 ) {
          for(int idx_cell = _pvtx_cell_idx[i_entity]; idx_cell < _pvtx_cell_idx[i_entity+1]; ++idx_cell) {
            int  i_cell = PDM_ABS(_pvtx_cell[idx_cell])-1;

            _pvtx_cell_coords[3*n_vtx_cell_to_send  ] = _pcell_center[3*i_cell  ];
            _pvtx_cell_coords[3*n_vtx_cell_to_send+1] = _pcell_center[3*i_cell+1];
            _pvtx_cell_coords[3*n_vtx_cell_to_send+2] = _pcell_center[3*i_cell+2];
            n_vtx_cell_to_send++;
          }
        }
      }
    }
    shift_part += fctv->n_part[i_domain];
  }

  PDM_distant_neighbor_exch(fctv->dn,
                            3 * sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            pvtx_cell_coords_n,
                  (void **) pvtx_cell_coords,
                           &fctv->pvtx_cell_coords_opp_n,
                 (void ***)&fctv->pvtx_cell_coords_opp);
  int    **pvtx_cell_coords_opp_n = fctv->pvtx_cell_coords_opp_n;
  double **pvtx_cell_coords_opp   = fctv->pvtx_cell_coords_opp;

  /* Same but for face_group */
  PDM_distant_neighbor_exch(fctv->dn,
                            3 * sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            fctv->vtx_face_bound_n,
                  (void **) fctv->vtx_face_bound_coords,
                           &fctv->pvtx_face_bound_coords_opp_n,
                 (void ***)&fctv->pvtx_face_bound_coords_opp);
  int    **pvtx_face_bound_coords_opp_n = fctv->pvtx_face_bound_coords_opp_n;
  double **pvtx_face_bound_coords_opp   = fctv->pvtx_face_bound_coords_opp;

  /*
   * Count receive
   */
  for(int i_part = 0; i_part < fctv->n_part_loc_all_domain; ++i_part){

    int nrecv = 0;
    int n_vtx = fctv->pn_vtx[i_part];
    int* _neighbor_idx       = fctv->neighbor_idx      [i_part];
    int* _neighbor_interface = fctv->neighbor_interface[i_part];

    for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
      for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {

        if(_neighbor_interface[idx_entity] != -40000) {
          int  i_interface = PDM_ABS(_neighbor_interface[idx_entity])-1;
          for(int idx_recv = 0; idx_recv < pvtx_cell_coords_opp_n[i_part][idx_entity]; ++idx_recv){
            for(int k = 0; k < 3; ++k) {
              pvtx_cell_coords_opp[i_part][3*(nrecv+idx_recv)+k] += PDM_SIGN(_neighbor_interface[idx_entity]) * fctv->pdi->translation_vect[i_interface][k];
            }
          }
        }
        nrecv += pvtx_cell_coords_opp_n[i_part][idx_entity];
      }
    }

    if(0 == 1) {
      printf("nrecv = %i \n",  nrecv);
      char filename[999];
      sprintf(filename, "opp_coords_%i_%i.vtk", i_rank, i_part);
      PDM_vtk_write_point_cloud(filename,
                                nrecv,
                                pvtx_cell_coords_opp[i_part],
                                NULL,
                                NULL);
    }

    /* Bnd */
    int nrecv_bnd = 0;

    for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
      for(int idx_entity = _neighbor_idx[i_entity]; idx_entity < _neighbor_idx[i_entity+1]; ++idx_entity) {

        if(_neighbor_interface[idx_entity] != -40000) {
          int  i_interface = PDM_ABS(_neighbor_interface[idx_entity])-1;
          for(int idx_recv = 0; idx_recv < pvtx_face_bound_coords_opp_n[i_part][idx_entity]; ++idx_recv){
            for(int k = 0; k < 3; ++k) {
              pvtx_face_bound_coords_opp[i_part][3*(nrecv_bnd+idx_recv)+k] += PDM_SIGN(_neighbor_interface[idx_entity]) * fctv->pdi->translation_vect[i_interface][k];
            }
          }
        }
        nrecv_bnd += pvtx_face_bound_coords_opp_n[i_part][idx_entity];
      }
    }

    if(0 == 1) {
      printf("nrecv_bnd = %i \n",  nrecv_bnd);
      char filename[999];
      sprintf(filename, "opp_coords_bnd_%i_%i.vtk", i_rank, i_part);
      PDM_vtk_write_point_cloud(filename,
                                nrecv_bnd,
                                pvtx_face_bound_coords_opp[i_part],
                                NULL,
                                NULL);
    }
  }



  /* Free */
  for(int i_part = 0; i_part < fctv->n_part_loc_all_domain; ++i_part){
    free(pvtx_cell_coords  [i_part]);
    free(pvtx_cell_coords_n[i_part]);
  }
  free(pvtx_cell_coords  );
  free(pvtx_cell_coords_n);

  /* Compute weight */



}



void
PDM_field_cell_to_vtx_part_set
(
  PDM_field_cell_to_vtx_t  *fctv,
  int                       i_domain,
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
  fctv->parts[i_domain][i_part].n_cell            = n_cell;
  fctv->parts[i_domain][i_part].n_face            = n_face;
  fctv->parts[i_domain][i_part].n_edge            = n_edge;
  fctv->parts[i_domain][i_part].n_vtx             = n_vtx;

  fctv->parts[i_domain][i_part].cell_face_idx = cell_face_idx;
  fctv->parts[i_domain][i_part].cell_face     = cell_face;

  fctv->parts[i_domain][i_part].face_edge_idx = face_edge_idx;
  fctv->parts[i_domain][i_part].face_edge     = face_edge;

  fctv->parts[i_domain][i_part].face_vtx_idx  = face_vtx_idx;
  fctv->parts[i_domain][i_part].face_vtx      = face_vtx;

  fctv->parts[i_domain][i_part].edge_vtx      = edge_vtx;

  fctv->parts[i_domain][i_part].cell_ln_to_gn = cell_ln_to_gn;
  fctv->parts[i_domain][i_part].face_ln_to_gn = face_ln_to_gn;
  fctv->parts[i_domain][i_part].edge_ln_to_gn = edge_ln_to_gn;
  fctv->parts[i_domain][i_part].vtx_ln_to_gn  = vtx_ln_to_gn;

  fctv->parts[i_domain][i_part].vtx = vtx_coord;
}

void
PDM_field_cell_to_vtx_part_mesh_nodal_set
(
  PDM_field_cell_to_vtx_t   *fctv,
  int                        i_domain,
  PDM_part_mesh_nodal_t     *pmn
)
{
  fctv->is_nodal = 1;
  fctv->pmn[i_domain] = pmn;
}

void
PDM_field_cell_to_vtx_graph_comm_set
(
  PDM_field_cell_to_vtx_t  *fctv,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  int                      *entity_part_bound_proc_idx,
  int                      *entity_part_bound_part_idx,
  int                      *entity_part_bound
)
{
  fctv->entity_part_bound_proc_idx[mesh_entity][i_domain][i_part] = entity_part_bound_proc_idx;
  fctv->entity_part_bound_part_idx[mesh_entity][i_domain][i_part] = entity_part_bound_part_idx;
  fctv->entity_part_bound         [mesh_entity][i_domain][i_part] = entity_part_bound;
  fctv->graph_comm_is_defined     [mesh_entity] = 1;
}


void
PDM_field_cell_to_vtx_part_domain_interface_shared_set
(
  PDM_field_cell_to_vtx_t     *fctv,
  PDM_part_domain_interface_t *pdi
)
{
  fctv->pdi = pdi;
}


void
PDM_field_cell_to_vtx_inverse_distance_weighting_p_set
(
  PDM_field_cell_to_vtx_t *fctv,
  int                      p
)
{
  fctv->idw_p = p;
}

void
PDM_field_cell_to_vtx_part_group_set
(
  PDM_field_cell_to_vtx_t  *fctv,
  int                       i_domain,
  int                       i_part,
  int                       i_group,
  PDM_bound_type_t          bound_type,
  int                       n_group_entity,
  int                      *group_entity
)
{
  assert(i_group < fctv->n_group[i_domain]);
  assert(bound_type == PDM_BOUND_TYPE_FACE);

  fctv->n_group_entity  [bound_type][i_domain][i_part][i_group] = n_group_entity;
  fctv->group_entity    [bound_type][i_domain][i_part][i_group] = group_entity;
  fctv->group_is_defined[bound_type] = 1;
}

void
PDM_field_cell_to_vtx_exch
(
        PDM_field_cell_to_vtx_t    *fctv,
        PDM_field_kind_t            field_kind,
        double                   ***local_field,
        double                  ****bound_field,
        double                  ****result_field
)
{
  int **pvtx_cell_n   = fctv->pvtx_cell_n;
  int **pvtx_cell_idx = fctv->pvtx_cell_idx;
  int **pvtx_cell     = fctv->pvtx_cell;

  int stride = _field_kind_to_size(field_kind);

  /*
   * Exchange volumic data
   */
  int    **pvtx_cell_field_n =  malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  double **pvtx_cell_field   =  malloc(fctv->n_part_loc_all_domain * sizeof(double *));
  int shift_part = 0;
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    /* First loop to count */
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
      int n_vtx          = fctv->pn_vtx      [i_part+shift_part];
      int* _neighbor_idx = fctv->neighbor_idx[i_part+shift_part];

      pvtx_cell_field_n[i_part+shift_part] = PDM_array_zeros_int(n_vtx);
      int *_pvtx_cell_field_n = pvtx_cell_field_n [i_part+shift_part];
      int *_pvtx_cell_n       = pvtx_cell_n       [i_part+shift_part];
      int *_pvtx_cell_idx     = pvtx_cell_idx     [i_part+shift_part];
      int *_pvtx_cell         = pvtx_cell         [i_part+shift_part];

      double*  _pcell_field = local_field[i_domain][i_part];

      /* Count */
      int n_vtx_cell_to_send = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        int n_neight = _neighbor_idx[i_entity+1] - _neighbor_idx[i_entity];
        if(n_neight > 0 ) {
          n_vtx_cell_to_send += _pvtx_cell_n[i_entity];
          _pvtx_cell_field_n[i_entity] = _pvtx_cell_n[i_entity];
        }
      }

      pvtx_cell_field[i_part+shift_part] = malloc( stride * n_vtx_cell_to_send * sizeof(double));
      double *_pvtx_cell_field = pvtx_cell_field[i_part+shift_part];

      n_vtx_cell_to_send = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        int n_neight = _neighbor_idx[i_entity+1] - _neighbor_idx[i_entity];
        if(n_neight > 0 ) {
          for(int idx_cell = _pvtx_cell_idx[i_entity]; idx_cell < _pvtx_cell_idx[i_entity+1]; ++idx_cell) {
            int  i_cell = PDM_ABS(_pvtx_cell[idx_cell])-1;

            for(int k = 0; k < stride; ++k){
              _pvtx_cell_field[stride*n_vtx_cell_to_send+k] = _pcell_field[stride*i_cell+k];
            }
            n_vtx_cell_to_send += 1; //stride;
          }
        }
      }
    }
    shift_part += fctv->n_part[i_domain];
  }

  /*
   * Exchange BND
   */
  shift_part = 0;
  int    **pvtx_face_field_n =  malloc(fctv->n_part_loc_all_domain * sizeof(int    *));
  double **pvtx_face_field   =  malloc(fctv->n_part_loc_all_domain * sizeof(double *));
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    /* First loop to count */
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
      int n_vtx          = fctv->pn_vtx      [i_part+shift_part];
      int* _neighbor_idx = fctv->neighbor_idx[i_part+shift_part];

      int *_vtx_face_bound_idx   = (int *) fctv->vtx_face_bound_idx  [i_part+shift_part];
      int *_vtx_face_bound_n     = (int *) fctv->vtx_face_bound_n    [i_part+shift_part];
      int *_vtx_face_bound       = (int *) fctv->vtx_face_bound      [i_part+shift_part];
      int *_vtx_face_bound_group = (int *) fctv->vtx_face_bound_group[i_part+shift_part];

      pvtx_face_field_n[i_part+shift_part] = PDM_array_zeros_int(n_vtx);
      int *_pvtx_face_field_n = pvtx_face_field_n [i_part+shift_part];

      /* Count */
      int n_vtx_face_to_send = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        int n_neight = _neighbor_idx[i_entity+1] - _neighbor_idx[i_entity];
        if(n_neight > 0 ) {
          n_vtx_face_to_send += _vtx_face_bound_n[i_entity];
          _pvtx_face_field_n[i_entity] = _vtx_face_bound_n[i_entity];
        }
      }

      pvtx_face_field[i_part+shift_part] = malloc( stride * n_vtx_face_to_send * sizeof(double));
      double  *_pvtx_face_field = pvtx_face_field[i_part+shift_part];
      double **_bound_field     = bound_field    [i_domain][i_part];

      n_vtx_face_to_send = 0;
      for(int i_entity = 0; i_entity < n_vtx; ++i_entity) {
        int n_neight = _neighbor_idx[i_entity+1] - _neighbor_idx[i_entity];
        if(n_neight > 0 ) {
          for(int idx_face = _vtx_face_bound_idx[i_entity]; idx_face < _vtx_face_bound_idx[i_entity+1]; ++idx_face) {
            int i_face  = _vtx_face_bound      [idx_face];
            int i_group = _vtx_face_bound_group[idx_face];

            for(int k = 0; k < stride; ++k){
              _pvtx_face_field[stride*n_vtx_face_to_send+k] = _bound_field[i_group][stride*i_face+k];
            }
            n_vtx_face_to_send += 1; //stride;
          }
        }
      }
    }
    shift_part += fctv->n_part[i_domain];
  }

  /*
   * Exchange
   */
  int    **pvtx_cell_field_opp_n = NULL;
  double **pvtx_cell_field_opp   = NULL;
  PDM_distant_neighbor_exch(fctv->dn,
                            stride * sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            pvtx_cell_field_n,
                  (void **) pvtx_cell_field,
                           &pvtx_cell_field_opp_n,
                 (void ***)&pvtx_cell_field_opp);

  /* Same but for face_group */
  int    **pvtx_face_field_opp_n = NULL;
  double **pvtx_face_field_opp   = NULL;
  PDM_distant_neighbor_exch(fctv->dn,
                            stride * sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            pvtx_face_field_n,
                  (void **) pvtx_face_field,
                           &pvtx_face_field_opp_n,
                 (void ***)&pvtx_face_field_opp);


  /*
   * Free  send
   */
  for(int i_part = 0; i_part < fctv->n_part_loc_all_domain; ++i_part ) {
    free(pvtx_cell_field_n[i_part]);
    free(pvtx_cell_field  [i_part]);
    free(pvtx_face_field_n[i_part]);
    free(pvtx_face_field  [i_part]);
  }
  free(pvtx_cell_field_n);
  free(pvtx_cell_field  );
  free(pvtx_face_field_n);
  free(pvtx_face_field  );

  /*
   * Post-treatment
   */
  _interpolate(fctv,
               field_kind,
               stride,
               local_field,
               bound_field,
               pvtx_cell_field_opp_n,
               pvtx_cell_field_opp,
               pvtx_face_field_opp_n,
               pvtx_face_field_opp,
               result_field);

  /*
   * Free
   */
  for(int i_part = 0; i_part < fctv->n_part_loc_all_domain; ++i_part ) {
    free(pvtx_cell_field_opp_n[i_part]);
    free(pvtx_cell_field_opp  [i_part]);
    free(pvtx_face_field_opp_n[i_part]);
    free(pvtx_face_field_opp  [i_part]);
  }
  free(pvtx_cell_field_opp_n);
  free(pvtx_cell_field_opp  );
  free(pvtx_face_field_opp_n);
  free(pvtx_face_field_opp  );


}

void
PDM_field_cell_to_vtx_set_cell_center
(
  PDM_field_cell_to_vtx_t  *fctv,
  int                       i_domain,
  int                       i_part,
  double                   *cell_center
)
{
  fctv->user_cell_center[i_domain][i_part] = cell_center;
}


void
PDM_field_cell_to_vtx_set_vtx_center
(
  PDM_field_cell_to_vtx_t  *fctv,
  int                       i_domain,
  int                       i_part,
  double                   *vtx_center
)
{
  fctv->user_vtx_center[i_domain][i_part] = vtx_center;
}


void
PDM_field_cell_to_vtx_free
(
 PDM_field_cell_to_vtx_t *fctv
)
{

  PDM_distant_neighbor_free(fctv->dn);

  for(int i_kind = 0; i_kind < PDM_MESH_ENTITY_MAX; ++i_kind) {
    fctv->graph_comm_is_defined     [i_kind] = 0;
    for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
      free(fctv->entity_part_bound_proc_idx[i_kind][i_domain]);
      free(fctv->entity_part_bound_part_idx[i_kind][i_domain]);
      free(fctv->entity_part_bound         [i_kind][i_domain]);
    }
    free(fctv->entity_part_bound_proc_idx[i_kind]);
    free(fctv->entity_part_bound_part_idx[i_kind]);
    free(fctv->entity_part_bound         [i_kind]);
  }
  free(fctv->entity_part_bound_proc_idx);
  free(fctv->entity_part_bound_part_idx);
  free(fctv->entity_part_bound         );
  free(fctv->graph_comm_is_defined     );

  for(int i_kind = 0; i_kind < PDM_GEOMETRY_KIND_MAX; ++i_kind) {
    fctv->group_is_defined     [i_kind] = 0;
    for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
        free(fctv->n_group_entity[i_kind][i_domain][i_part]);
        free(fctv->group_entity  [i_kind][i_domain][i_part]);
      }
      free(fctv->n_group_entity[i_kind][i_domain]);
      free(fctv->group_entity  [i_kind][i_domain]);
    }
    free(fctv->n_group_entity[i_kind]);
    free(fctv->group_entity  [i_kind]);
  }
  free(fctv->n_group_entity);
  free(fctv->group_entity);
  free(fctv->group_is_defined);


  int shift_part = 0;
  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < fctv->n_part[i_domain]; ++i_part) {
      free(fctv->neighbor_idx      [i_part+shift_part]);
      free(fctv->neighbor_desc     [i_part+shift_part]);
      free(fctv->neighbor_interface[i_part+shift_part]);

      free(fctv->vtx_face_bound_idx   [i_part+shift_part]);
      free(fctv->vtx_face_bound_n     [i_part+shift_part]);
      free(fctv->vtx_face_bound       [i_part+shift_part]);
      free(fctv->vtx_face_bound_group [i_part+shift_part]);
      free(fctv->vtx_face_bound_coords[i_part+shift_part]);

      free(fctv->pvtx_cell_n  [i_part+shift_part]);
      free(fctv->pvtx_cell_idx[i_part+shift_part]);
      free(fctv->pvtx_cell    [i_part+shift_part]);

      free(fctv->pvtx_cell_coords_opp_n      [i_part+shift_part]);
      free(fctv->pvtx_cell_coords_opp        [i_part+shift_part]);
      free(fctv->pvtx_face_bound_coords_opp_n[i_part+shift_part]);
      free(fctv->pvtx_face_bound_coords_opp  [i_part+shift_part]);
    }
    shift_part += fctv->n_part[i_domain];
  }
  free(fctv->neighbor_idx      );
  free(fctv->neighbor_desc     );
  free(fctv->neighbor_interface);
  free(fctv->pvtx_cell_n);
  free(fctv->pvtx_cell_idx);
  free(fctv->pvtx_cell);

  free(fctv->pvtx_cell_coords_opp_n      );
  free(fctv->pvtx_cell_coords_opp        );
  free(fctv->pvtx_face_bound_coords_opp_n);
  free(fctv->pvtx_face_bound_coords_opp  );

  free(fctv->vtx_face_bound_idx   );
  free(fctv->vtx_face_bound_n     );
  free(fctv->vtx_face_bound       );
  free(fctv->vtx_face_bound_group );
  free(fctv->vtx_face_bound_coords);
  free(fctv->pn_vtx);



  for(int i_domain = 0; i_domain < fctv->n_domain; ++i_domain) {
    free(fctv->parts      [i_domain]);
    free(fctv->cell_center[i_domain]);
  }
  free(fctv->parts);
  free(fctv->pmn);
  free(fctv->n_part_idx);
  free(fctv->n_part_g_idx);
  free(fctv->n_part);
  free(fctv->n_group);
  free(fctv->cell_center);

  free(fctv);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
