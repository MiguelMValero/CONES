
/*----------------------------------------------------------------------------
 *  Standar headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part_geom.h"
#include "pdm_hilbert.h"
#include "pdm_sort.h"
#include "pdm_geom_elem.h"
#include "pdm_binary_search.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_distrib.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dgeom_elem.h"
#include "pdm_array.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Compute cell center of elements
 *
 * \param [in]  comm         MPI Communicator
 * \param [in]  dn_cell       Number of cells in the current process
 * \param [in]  dn_face       Number of faces in the current process
 * \param [in]  dn_vtx        Number of vertices in the current process
 * \param [in]  dcell_face_idx Index of cell_face
 * \param [in]  dcell_face    cell face connectivity in the current process
 * \param [in]  dface_vtx_idx  Index of face_vtx
 * \param [in]  dface_vtx     face vertex connectivity in the current process
 * \param [in]  distrib_face    face distribution
 * \param [in]  dvtx_coord    coordinates of vertices
 * \param [in]  distrib_vtx     Vertex distribution
 *
 * \param [out] cell_center   Cell centers
 *
 */

void
PDM_dcompute_cell_center
(
  const PDM_MPI_Comm  comm,
  const int           dn_cell,
  const int          *dcell_face_idx,
  const PDM_g_num_t  *dcell_face,
  const int          *dface_vtx_idx,
  const PDM_g_num_t  *dface_vtx,
  const PDM_g_num_t  *distrib_face,
  const double       *dvtx_coord,
  const PDM_g_num_t  *distrib_vtx,
  double             *cell_center
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  /*PDM_log_trace_array_long (distrib_face, n_rank+1, "distrib_face : ");
    PDM_log_trace_array_long (distrib_vtx,  n_rank+1, "distrib_vtx  : ");*/

  int dn_face = (int) (distrib_face[i_rank+1] - distrib_face[i_rank]);

  PDM_g_num_t *dface_ln_to_gn = malloc (sizeof(PDM_g_num_t) * dn_face);
  for (int i = 0; i < dn_face; i++) {
    dface_ln_to_gn[i] = distrib_face[i_rank] + i + 1;
  }

  PDM_g_num_t *pvtx_ln_to_gn;
  int         *pface_vtx_idx;
  int         *pface_vtx;
  int          pn_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           dn_face,
                                                           dface_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pface_vtx_idx,
                                                           &pface_vtx);
  free (dface_ln_to_gn);

  /*PDM_log_trace_connectivity_long(dface_vtx_idx, dface_vtx, dn_face, "dface_vtx : ");
  PDM_log_trace_connectivity_int (pface_vtx_idx, pface_vtx, dn_face, "pface_vtx : ");
  for (int i = 0; i < pn_vtx; i++) {
    log_trace("vtx %d -> "PDM_FMT_G_NUM"\n", i+1, pvtx_ln_to_gn[i]);
  }*/

  double** tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_vtx,
                                        dvtx_coord,
                                        &pn_vtx,
                                        (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  double *pvtx_coord = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);
  free (pvtx_ln_to_gn);


  /* Compute face centers */
  double *dface_center = malloc (sizeof(double) * dn_face * 3);
  for (int i = 0; i < dn_face; i++) {
    for (int k = 0; k < 3; k++) {
      dface_center[3*i + k] = 0.;
    }

    double normalization = 1. / (double) (pface_vtx_idx[i+1] - pface_vtx_idx[i]);
    for (int j = pface_vtx_idx[i]; j < pface_vtx_idx[i+1]; j++) {
      int ivtx = pface_vtx[j] - 1;

      for (int k = 0; k < 3; k++) {
        dface_center[3*i + k] += pvtx_coord[3*ivtx + k];
      }
    }

    for (int k = 0; k < 3; k++) {
      dface_center[3*i + k] *= normalization;
    }
  }
  free (pvtx_coord);
  free (pface_vtx_idx);
  free (pface_vtx);

  /* Compute cell centers */
  PDM_compute_center_from_descending_connectivity (dcell_face_idx,
                                                   dcell_face,
                                                   dn_cell,
                                                   distrib_face,
                                                   cell_center,
                                                   dface_center,
                                                   comm);
  free (dface_center);
}
/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_part_entity_geom
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const PDM_g_num_t   dn_entity,
 const double       *dentity_coord,
 const double       *dentity_weight,
       int          *dentity_part
)
{
  PDM_UNUSED(method);
  assert (method == PDM_PART_GEOM_HILBERT);

  const int dim = 3;

  /** TRAITEMENT HILBERT FVM **/
  PDM_hilbert_code_t *hilbert_codes     = (PDM_hilbert_code_t *) malloc (dn_entity * sizeof(PDM_hilbert_code_t));

  /** Initialisation **/

  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, dn_entity, dentity_coord, extents, comm);

  /** Hilbert Coordinates Computation **/
  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, dn_entity, dentity_coord, hilbert_codes);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  int n_total_part;
  PDM_MPI_Allreduce ((void *) &n_part, &n_total_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  PDM_hilbert_code_t *hilbert_codes_idx = (PDM_hilbert_code_t *) malloc ((n_total_part+1) * sizeof(PDM_hilbert_code_t));

  double * weight = (double *) malloc (dn_entity * sizeof(double));
  if (dentity_weight != NULL) {
    for(int i = 0; i < dn_entity; ++i) {
      weight [i] = dentity_weight [i];
    }
  }
  else {
    for(int i = 0; i < dn_entity; ++i) {
      weight [i] = 1.;
    }
  }

  PDM_hilbert_build_rank_index (dim,
                                n_total_part,
                                dn_entity,
                                hilbert_codes,
                                weight,
                                NULL,
                                hilbert_codes_idx,
                                comm);

  free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/

  for(int i = 0; i < dn_entity; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_total_part,
                                                hilbert_codes[i],
                                                hilbert_codes_idx);
    dentity_part[i] = (int) quantile;

  }

  free(hilbert_codes_idx);
  free(hilbert_codes);
}

/**
 *
 * \brief Perform geomtric patitionning
 *
 * \param [in]   method         Geometric method
 * \param [in]   n_part          Number of partition to build on this process
 * \param [in]   comm           Communicator
 * \param [in]   dn_cell         Number of distributed cells
 * \param [in]   dn_face         Number of distributed faces
 * \param [in]   dn_vtx          Number of distributed vertices
 * \param [in]   dcell_face_idx   Distributed cell face connectivity index or NULL
 *                              (size : dn_cell + 1, numbering : 0 to n-1)
 * \param [in]   dcell_face      Distributed cell face connectivity or NULL
 *                              (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
 * \param [in]   dcell_weight    Cell weight (size : n_cell) or NULL
 * \param [in]   dface_vtx_idx    Distributed face to vertex connectivity index
 *                              (size : dn_face + 1, numbering : 0 to n-1)
 * \param [in]   dface_vtx       Distributed face to vertex connectivity
 *                              (size : dface_vtx_idx[dn_face], numbering : 1 to n)
 * \param [in]   dvtx_coord      Distributed vertex coordinates
 *                              (size : 3*dn_vtx)
 * \param [inout]   dcell_part      Distributed cell partitioning
 *                              (size = dn_cell)
 *
 */

void
PDM_part_geom
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_cell,
 const int          *dcell_face_idx,
 const PDM_g_num_t  *dcell_face,
 const int          *dcell_weight,
 const int          *dface_vtx_idx,
 const PDM_g_num_t  *dface_vtx,
 const PDM_g_num_t  *distrib_face,
 const double       *dvtx_coord,
 const PDM_g_num_t  *distrib_vtx,
       int          *dcell_part
)
{
  assert (method == PDM_PART_GEOM_HILBERT);
  /*
   * cell center computation
   */
  double *barycenter_coords = (double *) malloc (dn_cell * 3 * sizeof(double ));
  PDM_dcompute_cell_center (comm,
                            dn_cell,
                            dcell_face_idx,
                            dcell_face,
                            dface_vtx_idx,
                            dface_vtx,
                            distrib_face,
                            dvtx_coord,
                            distrib_vtx,
                            barycenter_coords);

  double *dcell_weight_d = NULL;
  if(dcell_weight != NULL) {
    dcell_weight_d = malloc(dn_cell * sizeof(double));
    for(int i = 0; i < dn_cell; ++i) {
      dcell_weight_d[i] = dcell_weight[i];
    }
  }

  PDM_part_entity_geom(method,
                       n_part,
                       comm,
                       dn_cell,
                       barycenter_coords,
                       dcell_weight_d,
                       dcell_part);

  if(dcell_weight != NULL) {
    free(dcell_weight_d);
  }

  free(barycenter_coords);
}


void
PDM_part_geom_0d
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_vtx,
 const double       *dvtx_coord,
 const double       *dvtx_weight,
       int          *dvtx_part
)
{
  PDM_part_entity_geom(method,
                       n_part,
                       comm,
                       dn_vtx,
                       dvtx_coord,
                       dvtx_weight,
                       dvtx_part);
}

void
PDM_part_geom_1d
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_edge,
 const int           dn_vtx,
 const PDM_g_num_t  *dedge_vtx,
 const double       *dvtx_coord,
 const double       *dedge_weight,
       int          *dedge_part
)
{
  PDM_g_num_t *distrib_vtx = PDM_compute_entity_distribution(comm, dn_vtx);

  int *dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);

  double *dedge_center = (double *) malloc(sizeof(double) * dn_edge * 3);
  PDM_compute_center_from_descending_connectivity(dedge_vtx_idx,
                                                  dedge_vtx,
                                                  dn_edge,
                                                  distrib_vtx,
                                                  dedge_center,
                                  (double *)      dvtx_coord,
                                                  comm);

  PDM_part_entity_geom(method,
                       n_part,
                       comm,
                       dn_edge,
                       dedge_center,
                       dedge_weight,
                       dedge_part);

  free(distrib_vtx);
  free(dedge_center);
  free(dedge_vtx_idx);
}

void
PDM_part_geom_2d
(
 PDM_part_geom_t     method,
 const int           n_part,
 const PDM_MPI_Comm  comm,
 const int           dn_face,
 const int           dn_edge,
 const int           dn_vtx,
 const int          *dface_vtx_idx,
 const PDM_g_num_t  *dface_vtx,
 const int          *dface_edge_idx,
 const PDM_g_num_t  *dface_edge,
 const PDM_g_num_t  *dedge_vtx,
 const double       *dvtx_coord,
 const double       *dface_weight,
       int          *dface_part
)
{

  PDM_g_num_t *distrib_vtx = PDM_compute_entity_distribution(comm, dn_vtx);
  double *dface_center = (double *) malloc(sizeof(double) * dn_face * 3);

  if(dface_vtx_idx != NULL) {
    PDM_compute_center_from_descending_connectivity(dface_vtx_idx,
                                                    dface_vtx,
                                                    dn_face,
                                                    distrib_vtx,
                                                    dface_center,
                                    (double *)      dvtx_coord,
                                                    comm);

  } else {
    assert(dface_edge_idx != NULL);
    int *dedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, dn_edge);

    double *dedge_center = (double *) malloc(sizeof(double) * dn_edge * 3);

    PDM_compute_center_from_descending_connectivity(dedge_vtx_idx,
                                                    dedge_vtx,
                                                    dn_edge,
                                                    distrib_vtx,
                                                    dedge_center,
                                    (double *)      dvtx_coord,
                                                    comm);

    PDM_g_num_t *distrib_edge = PDM_compute_entity_distribution(comm, dn_edge);
    PDM_compute_center_from_descending_connectivity(dface_edge_idx,
                                                    dface_edge,
                                                    dn_face,
                                                    distrib_edge,
                                                    dface_center,
                                    (double *)      dedge_center,
                                                    comm);


    free(dedge_vtx_idx);
    free(dedge_center);
    free(distrib_edge);
  }

  PDM_part_entity_geom(method,
                       n_part,
                       comm,
                       dn_face,
                       dface_center,
                       dface_weight,
                       dface_part);

  free(distrib_vtx);
  free(dface_center);
}


void
PDM_dreorder_from_coords
(
 PDM_part_geom_t  method,
 int              dim,
 PDM_g_num_t     *distrib_vtx,
 double          *dcoords,
 PDM_g_num_t     *ln_to_gn,
 PDM_MPI_Comm     comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  assert (method == PDM_PART_GEOM_HILBERT);

  int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];

  PDM_hilbert_code_t *hilbert_codes     = (PDM_hilbert_code_t *) malloc (dn_vtx * sizeof(PDM_hilbert_code_t));

  /** Initialisation **/
  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, dn_vtx, dcoords, extents, comm);

  /** Hilbert Coordinates Computation **/
  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, dn_vtx, dcoords, hilbert_codes);

  // PDM_hilbert_code_t *tmp_hilbert_codes = (PDM_hilbert_code_t *) malloc (dn_vtx * sizeof(PDM_hilbert_code_t));
  // for (int i = 0; i < dn_vtx; ++i) {
  //   tmp_hilbert_codes [i] = hilbert_codes [i];
  // }

  ///** Calcul des index des codes Hilbert **/
  // int * hilbert_order = (int * ) malloc (dn_vtx * sizeof(int));
  // for (int i = 0; i < dn_vtx; ++i) {
  //   hilbert_order [i] = i;
  // }
  // assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  // PDM_sort_double (tmp_hilbert_codes, NULL, dn_vtx);
  // free(tmp_hilbert_codes);

  PDM_hilbert_code_t *hilbert_codes_idx = (PDM_hilbert_code_t *) malloc ((n_rank+1) * sizeof(PDM_hilbert_code_t));

  double *weight = PDM_array_const_double(dn_vtx, 1.);
  PDM_hilbert_build_rank_index (dim,
                                n_rank,
                                dn_vtx,
                                hilbert_codes,
                                weight,
                                NULL, // No need order
                                hilbert_codes_idx,
                                comm);
  free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/
  for(int i = 0; i < dn_vtx; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_rank,
                                                  hilbert_codes[i],
                                                  hilbert_codes_idx);
    ln_to_gn [i] = (int) (quantile + 1); // Because ln_to_gn of part_to_block begin at 1
  }

  // part_to_block avec ln_to_gn 1 2 3 4 .... pdm_assembly_partition
  // Puis on échange les hilbert_codes, qu'on retrie localement

  PDM_g_num_t* distrib_rank = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_rank+1; ++i) {
    distrib_rank[i] = i;
  }

  /*
   * Each proc get all the entities affected to its partitions
   */
  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_NOTHING,
                             1.,
                             &ln_to_gn,
                             distrib_rank,
                             &dn_vtx,
                             1,
                             comm);
  free(distrib_rank);

  const int n_vtx_block = PDM_part_to_block_n_elt_block_get (ptb);

  // log_trace("n_vtx_block = %i | dn_vtx = %i \n", n_vtx_block, dn_vtx);

  double *blk_hilbert_codes = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &hilbert_codes,
                          NULL,
                (void **) &blk_hilbert_codes);

  /* Resend */
  for(int i = 0; i < dn_vtx; ++i) {
    ln_to_gn [i] = distrib_vtx[i_rank] + i + 1; // Donc correspond a la numeration absolu initiale
  }

  PDM_g_num_t* blk_ln_to_gn;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &ln_to_gn,
                          NULL,
                (void **) &blk_ln_to_gn);

  free(hilbert_codes_idx);
  free(hilbert_codes);
  PDM_part_to_block_free(ptb);

  /* Reorder locally */
  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  int *hilbert_order = (int * ) malloc( n_vtx_block * sizeof(int));
  for (int i = 0; i < n_vtx_block; ++i) {
    hilbert_order [i] = i;
  }
  PDM_sort_double (blk_hilbert_codes, hilbert_order, n_vtx_block);
  free(blk_hilbert_codes);


  /* Apply order to blk_ln_to_gn */
  PDM_g_num_t* sorted_blk_ln_to_gn = malloc( n_vtx_block * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_vtx_block; ++i) {
    sorted_blk_ln_to_gn[i] = blk_ln_to_gn[hilbert_order[i]];
  }
  free(blk_ln_to_gn);
  free(hilbert_order);

  PDM_g_num_t* distrib_blk_vtx = PDM_compute_entity_distribution(comm, n_vtx_block);
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_blk_vtx,
                              (const PDM_g_num_t **)  &ln_to_gn,
                                                      &dn_vtx,
                                                      1,
                                                      comm);

  int stride_one = 1;
  PDM_block_to_part_exch_in_place(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
              (void *)   sorted_blk_ln_to_gn,
                         NULL,
              (void **) &ln_to_gn);
  PDM_block_to_part_free(btp);
  free(sorted_blk_ln_to_gn);
  free(distrib_blk_vtx);
}


void
PDM_dreorder_from_length
(
 int              dim,
 PDM_g_num_t     *distrib_in,
 double          *length,
 PDM_g_num_t     *ln_to_gn,
 PDM_MPI_Comm     comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dn_length = distrib_in[i_rank+1] - distrib_in[i_rank];

  // PDM_hilbert_code_t *hilbert_codes     = (PDM_hilbert_code_t *) malloc (dn_length * sizeof(PDM_hilbert_code_t));
  PDM_hilbert_code_t *tmp_hilbert_codes = (PDM_hilbert_code_t *) malloc (dn_length * sizeof(PDM_hilbert_code_t));

  for (int i = 0; i < dn_length; ++i) {
    tmp_hilbert_codes [i] = length [i];
  }

  ///** Calcul des index des codes Hilbert **/
  int * hilbert_order = (int * ) malloc (dn_length * sizeof(int));

  for (int i = 0; i < dn_length; ++i) {
    hilbert_order [i] = i;
  }

  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  PDM_sort_double (tmp_hilbert_codes, hilbert_order, dn_length);

  free(tmp_hilbert_codes);

  PDM_hilbert_code_t *hilbert_codes_idx = (PDM_hilbert_code_t *) malloc ((n_rank+1) * sizeof(PDM_hilbert_code_t));

  double *weight = PDM_array_const_double(dn_length, 1.);
  PDM_hilbert_build_rank_index (dim,
                                n_rank,
                                dn_length,
                                length,
                                weight,
                                hilbert_order,
                                hilbert_codes_idx,
                                comm);


  free(weight);

  /** Remplissage de cell_parts -> en fct des codes Hilbert **/
  for(int i = 0; i < dn_length; ++i) {
    size_t quantile = PDM_hilbert_quantile_search(n_rank,
                                                  length[i],
                                                  hilbert_codes_idx);
    ln_to_gn [i] = (int) (quantile + 1); // Because ln_to_gn of part_to_block begin at 1
  }

  // part_to_block avec ln_to_gn 1 2 3 4 .... pdm_assembly_partition
  // Puis on échange les hilbert_codes, qu'on retrie localement

  PDM_g_num_t* distrib_rank = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_rank+1; ++i) {
    distrib_rank[i] = i;
  }

  /*
   * Each proc get all the entities affected to its partitions
   */
  PDM_part_to_block_t *ptb =
   PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                             PDM_PART_TO_BLOCK_POST_NOTHING,
                             1.,
                             &ln_to_gn,
                             distrib_rank,
                             &dn_length,
                             1,
                             comm);
  free(distrib_rank);

  const int n_vtx_block = PDM_part_to_block_n_elt_block_get (ptb);

  // log_trace("n_vtx_block = %i | dn_length = %i \n", n_vtx_block, dn_length);

  double *blk_hilbert_codes = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &length,
                          NULL,
                (void **) &blk_hilbert_codes);

  /* Resend */
  for(int i = 0; i < dn_length; ++i) {
    ln_to_gn [i] = distrib_in[i_rank] + i + 1; // Donc correspond a la numeration absolu initiale
  }

  PDM_g_num_t* blk_ln_to_gn;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void **) &ln_to_gn,
                          NULL,
                (void **) &blk_ln_to_gn);

  free(hilbert_codes_idx);
  PDM_part_to_block_free(ptb);

  /* Reorder locally */
  assert (sizeof(double) == sizeof(PDM_hilbert_code_t));
  hilbert_order = (int * ) realloc(hilbert_order,  n_vtx_block * sizeof(int));
  for (int i = 0; i < n_vtx_block; ++i) {
    hilbert_order [i] = i;
  }
  PDM_sort_double (blk_hilbert_codes, hilbert_order, n_vtx_block);
  //PDM_log_trace_array_double(blk_hilbert_codes, n_vtx_block, "tmp_edge_length : ");
  free(blk_hilbert_codes);


  /* Apply order to blk_ln_to_gn */
  PDM_g_num_t* sorted_blk_ln_to_gn = malloc( n_vtx_block * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_vtx_block; ++i) {
    sorted_blk_ln_to_gn[i] = blk_ln_to_gn[hilbert_order[i]];
  }
  free(blk_ln_to_gn);
  free(hilbert_order);

  PDM_g_num_t* distrib_blk_vtx = PDM_compute_entity_distribution(comm, n_vtx_block);
  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_blk_vtx,
                              (const PDM_g_num_t **)  &ln_to_gn,
                                                      &dn_length,
                                                      1,
                                                      comm);

  int stride_one = 1;
  PDM_block_to_part_exch_in_place(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
              (void *)   sorted_blk_ln_to_gn,
                         NULL,
              (void **) &ln_to_gn);
  PDM_block_to_part_free(btp);
  free(sorted_blk_ln_to_gn);
  free(distrib_blk_vtx);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
