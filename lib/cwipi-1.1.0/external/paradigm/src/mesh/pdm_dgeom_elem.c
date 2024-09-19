/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_polygon.h"
#include "pdm_plane.h"
#include "pdm_dgeom_elem.h"
#include "pdm_block_to_part.h"
#include "pdm_plane.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_dconnectivity_transform.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_compute_center_from_descending_connectivity
(
  const int         *dentity1_entity2_idx,
  const PDM_g_num_t *dentity1_entity2,
  const int          dn_entity1,
  const PDM_g_num_t *dentity2_distrib,
  double            *dentity1_coord,
  double            *dentity2_coord,
  PDM_MPI_Comm       comm
)
{
  int *dentity1_entity2_sgn = malloc(dentity1_entity2_idx[dn_entity1] * sizeof(int));
  PDM_g_num_t *dentity1_entity2_abs = malloc(dentity1_entity2_idx[dn_entity1] * sizeof(PDM_g_num_t));
  for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
    dentity1_entity2_sgn[i] = PDM_SIGN(dentity1_entity2[i]);
    dentity1_entity2_abs[i] = PDM_ABS (dentity1_entity2[i]);
  }
  PDM_block_to_part_t *btp_entity1_coord = PDM_block_to_part_create (dentity2_distrib,
                                              (const PDM_g_num_t **) &dentity1_entity2_abs,
                                                                     &dentity1_entity2_idx[dn_entity1],
                                                                     1,
                                                                     comm);
  /*for(int i = 0; i < dentity1_entity2_idx[dn_entity1]; ++i) {
    dentity1_entity2[i]     = dentity1_entity2[i]*dentity1_entity2_sgn[i];
    }*/
  free(dentity1_entity2_sgn);
  free(dentity1_entity2_abs);

  int strid_one = 1;
  double **tmp_entity1_entity2_coord;
  PDM_block_to_part_exch (btp_entity1_coord,
                           3 * sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           &strid_one,
                  (void *) dentity2_coord,
                           NULL,
                (void ***) &tmp_entity1_entity2_coord);
  double *dentity1_entity2_coord = tmp_entity1_entity2_coord[0];
  free(tmp_entity1_entity2_coord);
  PDM_block_to_part_free(btp_entity1_coord);

  for(int i_entity1 = 0; i_entity1 < dn_entity1; ++i_entity1) {
    dentity1_coord[3*i_entity1  ] = 0.;
    dentity1_coord[3*i_entity1+1] = 0.;
    dentity1_coord[3*i_entity1+2] = 0.;
    int n_entity2_per_entity1 = dentity1_entity2_idx[i_entity1+1] - dentity1_entity2_idx[i_entity1];
    double inv = 1./n_entity2_per_entity1;
    for(int idx_entity2 = dentity1_entity2_idx[i_entity1]; idx_entity2 < dentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {
      dentity1_coord[3*i_entity1  ] += dentity1_entity2_coord[3*idx_entity2  ];
      dentity1_coord[3*i_entity1+1] += dentity1_entity2_coord[3*idx_entity2+1];
      dentity1_coord[3*i_entity1+2] += dentity1_entity2_coord[3*idx_entity2+2];
    }
    dentity1_coord[3*i_entity1  ] = dentity1_coord[3*i_entity1  ]*inv;
    dentity1_coord[3*i_entity1+1] = dentity1_coord[3*i_entity1+1]*inv;
    dentity1_coord[3*i_entity1+2] = dentity1_coord[3*i_entity1+2]*inv;

  }
  free(dentity1_entity2_coord);

}


void
PDM_compute_dface_normal
(
  const int         *dface_vtx_idx,
  const PDM_g_num_t *dface_vtx,
  const int          dn_face,
  const PDM_g_num_t *dvtx_distrib,
  double            *dvtx_coord,
  double            *dface_normal,
  PDM_MPI_Comm       comm
)
{
  int *dface_vtx_sgn = malloc(dface_vtx_idx[dn_face] * sizeof(int));
  PDM_g_num_t *dface_vtx_abs = malloc(dface_vtx_idx[dn_face] * sizeof(PDM_g_num_t));
  for(int i = 0; i < dface_vtx_idx[dn_face]; ++i) {
    dface_vtx_sgn[i] = PDM_SIGN(dface_vtx[i]);
    dface_vtx_abs[i] = PDM_ABS (dface_vtx[i]);
  }
  PDM_block_to_part_t *btp_entity1_coord = PDM_block_to_part_create (dvtx_distrib,
                                              (const PDM_g_num_t **) &dface_vtx_abs,
                                                                     &dface_vtx_idx[dn_face],
                                                                     1,
                                                                     comm);
  free(dface_vtx_sgn);
  free(dface_vtx_abs);

  int strid_one = 1;
  double **tmp_face_vtx_coord;
  PDM_block_to_part_exch (btp_entity1_coord,
                           3 * sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           &strid_one,
                  (void *) dvtx_coord,
                           NULL,
                (void ***) &tmp_face_vtx_coord);
  double *dface_vtx_coord = tmp_face_vtx_coord[0];
  free(tmp_face_vtx_coord);
  PDM_block_to_part_free(btp_entity1_coord);

  // double* _dface_vtx_ptr = dface_vtx_coord;
  for(int i_face = 0; i_face < dn_face; ++i_face) {

    dface_normal[3*i_face  ] = i_face;
    dface_normal[3*i_face+1] = i_face;
    dface_normal[3*i_face+2] = i_face;

    int n_vtx_per_face = dface_vtx_idx[i_face+1] - dface_vtx_idx[i_face];
    PDM_plane_normal(n_vtx_per_face, dface_vtx_coord + 3*dface_vtx_idx[i_face], &dface_normal[3*i_face  ]);

    // _dface_vtx_ptr += 3 * n_vtx_per_face;

  }
  free(dface_vtx_coord);
}


void
PDM_compute_vtx_characteristic_length
(
 PDM_MPI_Comm    comm,
 int             dn_face,
 int             dn_edge,
 int             dn_vtx,
 int            *dface_vtx_idx,
 PDM_g_num_t    *dface_vtx,
 PDM_g_num_t    *dedge_vtx,
 double         *dvtx_coord,
 double        **dchar_length_out
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *distrib_vtx  = PDM_compute_entity_distribution(comm, dn_vtx );
  PDM_g_num_t *distrib_edge = PDM_compute_entity_distribution(comm, dn_edge);
  PDM_g_num_t *distrib_face = PDM_compute_entity_distribution(comm, dn_face);

  // Compute graph of vtx
  PDM_g_num_t *dvtx_vtx_idx = NULL;
  PDM_g_num_t *dvtx_vtx     = NULL;
  if(dedge_vtx == NULL) {
    int         *dvtx_face_idx = NULL;
    PDM_g_num_t *dvtx_face     = NULL;
    PDM_dconnectivity_transpose(comm,
                                distrib_face,
                                distrib_vtx,
                                dface_vtx_idx,
                                dface_vtx,
                                0,
                                &dvtx_face_idx,
                                &dvtx_face);

    PDM_deduce_combine_connectivity_dual(comm,
                                         distrib_vtx,
                                         distrib_face,
                                         dvtx_face_idx,
                                         dvtx_face,
                                         dface_vtx_idx,
                                         dface_vtx,
                                         0,
                                         &dvtx_vtx_idx,
                                         &dvtx_vtx);
    free(dvtx_face_idx);
    free(dvtx_face    );

  } else {

    int *dedge_vtx_idx = malloc((dn_edge+1) * sizeof(int));
    for(int i = 0; i < dn_edge+1; ++i) {
      dedge_vtx_idx[i] = 2 * i;
    }
    int         *dvtx_edge_idx = NULL;
    PDM_g_num_t *dvtx_edge     = NULL;
    PDM_dconnectivity_transpose(comm,
                                distrib_edge,
                                distrib_vtx,
                                dedge_vtx_idx,
                                dedge_vtx,
                                0,
                                &dvtx_edge_idx,
                                &dvtx_edge);

    PDM_deduce_combine_connectivity_dual(comm,
                                         distrib_vtx,
                                         distrib_edge,
                                         dvtx_edge_idx,
                                         dvtx_edge,
                                         dedge_vtx_idx,
                                         dedge_vtx,
                                         0,
                                         &dvtx_vtx_idx,
                                         &dvtx_vtx);

    free(dedge_vtx_idx);
    free(dvtx_edge_idx);
    free(dvtx_edge    );
  }

  /*
   * Partitionnement du pauvre
   */
  PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_vtx,
                              (const PDM_g_num_t **)  &dvtx_vtx,
                              (const int *)           &dvtx_vtx_idx[dn_vtx],
                                                      1,
                                                      comm);

  int stride_one = 1;
  double **tmp_vtx_vtx_coord = NULL;
  PDM_block_to_part_exch(btp,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         dvtx_coord,
                         NULL,
          (void ***)    &tmp_vtx_vtx_coord);
  double *pvtx_vtx_coord = tmp_vtx_vtx_coord[0];
  free(tmp_vtx_vtx_coord);
  PDM_block_to_part_free(btp);
  free(dvtx_vtx    );

  double *char_length = malloc(dn_vtx * sizeof(double));

  for(int i_vtx = 0; i_vtx < dn_vtx; ++i_vtx) {

    char_length[i_vtx] = HUGE_VAL;

    for(int idx_vtx = dvtx_vtx_idx[i_vtx]; idx_vtx < dvtx_vtx_idx[i_vtx+1]; ++idx_vtx) {

      double length2 = 0.;
      for(int k = 0; k < 3; ++k) {
        double delta = dvtx_coord[3*i_vtx + k] - pvtx_vtx_coord[3*idx_vtx + k];
        length2 += delta * delta;
      }

      char_length[i_vtx] = PDM_MIN(char_length[i_vtx], length2);
    }

    // char_length[i_vtx] = PDM_MAX(eps_base, tol*sqrt(char_length[i_vtx]));
    char_length[i_vtx] = sqrt(char_length[i_vtx]);
  }

  *dchar_length_out = char_length;

  free(distrib_vtx );
  free(distrib_edge);
  free(distrib_face);
  free(pvtx_vtx_coord);
  free(dvtx_vtx_idx);
}





#ifdef __cplusplus
}
#endif /* __cplusplus */
