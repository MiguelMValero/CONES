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
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_array.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_multipart.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
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
           int          *elt_type,
           int          *n_steps)
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
    else if (strcmp(argv[i], "-n_steps") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_steps = atoi(argv[i]);
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static inline double _rand (void)
{
  return 2 * (double) rand() / (double) RAND_MAX - 1;
}


static void
_generate_mesh
(
 PDM_MPI_Comm          comm,
const PDM_g_num_t      n_vtx_seg,
PDM_Mesh_nodal_elt_t   elt_type,
PDM_multipart_t      **_mpart
 )
{
  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 2);

  double length = 1.;
  int n_part = 1;
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         0,
                                                         0,
                                                         0.,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
  double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
  int dn_vtx = (int) (vtx_distrib[i_rank+1] - vtx_distrib[i_rank]);

  for (int i = 0; i < 2*vtx_distrib[i_rank]; i++) {
    rand();
  }

  double noise = 0.49 * length / (double) (n_vtx_seg - 1);
  for (int i = 0; i < dn_vtx; i++) {
    double x = dvtx_coord[3*i];
    double y = dvtx_coord[3*i+1];

    for (int j = 0; j < 2; j++) {
      double dx = noise * _rand();
      if (x > 0. && x < length && y > 0 && y < length) {
        dvtx_coord[3*i+j] += dx;
      }
    }
  }



  int n_domain = 1;
  int *n_part_domains = (int *) malloc(sizeof(int) * n_domain);
  n_part_domains[0] = n_part;

  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                n_part_domains,
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

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);


  /* Run */
  PDM_multipart_compute (mpart);

  free(n_part_domains);

  *_mpart = mpart;

  PDM_dcube_nodal_gen_free(dcube);
}




// static void
// _get_groups
// (
//  PDM_multipart_t     *multipart,
//  const int            i_domain,
//  const int            i_part,
//        int           *n_face_group,
//        int          **face_bound_idx,
//        int          **face_bound
//  )
// {

//   int n_section;
//   int*n_elt;
//   int n_cell;
//   int n_face;
//   int n_face_part_bound;
//   int n_vtx;
//   int n_proc;
//   int n_total_part;
//   int s_cell_face;
//   int s_face_vtx;
//   int s_face_bound;
//   int s_face_join;
//   int n_join_groups;
//   PDM_multipart_part_dim_get(multipart,
//                              i_domain,
//                              i_part,
//                              &n_section,
//                              &n_elt,
//                              &n_cell,
//                              &n_face,
//                              &n_face_part_bound,
//                              &n_vtx,
//                              &n_proc,
//                              &n_total_part,
//                              &s_cell_face,
//                              &s_face_vtx,
//                              &s_face_bound,
//                              n_face_group,
//                              &s_face_join,
//                              &n_join_groups);


//   int         **elt_vtx_idx;
//   int         **elt_vtx;
//   PDM_g_num_t **elt_section_ln_to_gn;
//   int          *cell_tag;
//   int          *cell_face_idx;
//   int          *cell_face;
//   PDM_g_num_t  *cell_ln_to_gn;
//   int          *face_tag;
//   int          *face_cell;
//   int          *face_vtx_idx;
//   int          *face_vtx;
//   PDM_g_num_t  *face_ln_to_gn;
//   int          *face_part_bound_proc_idx;
//   int          *face_part_bound_part_idx;
//   int          *face_part_bound;
//   int          *vtx_tag;
//   double       *vtx;
//   PDM_g_num_t  *vtx_ln_to_gn;
//   PDM_g_num_t  *face_bound_ln_to_gn;
//   int          *face_join_idx;
//   int          *face_join;
//   PDM_g_num_t  *face_join_ln_to_gn;
//   PDM_multipart_part_val_get(multipart,
//                              i_domain,
//                              i_part,
//                              &elt_vtx_idx,
//                              &elt_vtx,
//                              &elt_section_ln_to_gn,
//                              &cell_tag,
//                              &cell_face_idx,
//                              &cell_face,
//                              &cell_ln_to_gn,
//                              &face_tag,
//                              &face_cell,
//                              &face_vtx_idx,
//                              &face_vtx,
//                              &face_ln_to_gn,
//                              &face_part_bound_proc_idx,
//                              &face_part_bound_part_idx,
//                              &face_part_bound,
//                              &vtx_tag,
//                              &vtx,
//                              &vtx_ln_to_gn,
//                              face_bound_idx,
//                              face_bound,
//                              &face_bound_ln_to_gn,
//                              &face_join_idx,
//                              &face_join,
//                              &face_join_ln_to_gn);
// }




/**
 *
 * \brief  Main
 *
 */
int main(int argc, char *argv[])
{

  /*
   *
   * elt_type :
   *  2 -> tria
   *  3 -> quad
   */

  PDM_g_num_t          n_vtx_seg = 5;  // Number of vtx on each side of the cube mesh
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TRIA3;  // Type of cells
  int                  n_steps   = 10; // Number of smoothing steps
  _read_args(argc,
             argv,
             &n_vtx_seg,
     (int *) &elt_type,
             &n_steps);


  PDM_MPI_Init (&argc, &argv);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_multipart_t *mpart = NULL;
  _generate_mesh(comm,
                 n_vtx_seg,
                 elt_type,
                 &mpart);

  // int          pn_face         = 0;
  int          pn_edge         = 0;
  int          pn_vtx          = 0;
  int         *pface_edge_idx  = NULL;
  int         *pface_edge      = NULL;
  int         *pedge_vtx       = NULL;
  double      *pvtx_coord      = NULL;
  int          pn_edge_group   = 0;
  int         *pgroup_edge_idx = NULL;
  int         *pgroup_edge     = NULL;
  PDM_g_num_t *vtx_ln_to_gn   = NULL;

  /* Get vertices */
  int i_part = 0;
  pn_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                           0,
                                           i_part,
                                           &pvtx_coord,
                                           PDM_OWNERSHIP_KEEP);

  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_VTX,
                                  &vtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  /* Get edges */
  int* tmp_pedge_vtx_idx = NULL;
  pn_edge = PDM_multipart_part_connectivity_get(mpart,
                                               0,
                                               i_part,
                                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                               &tmp_pedge_vtx_idx,
                                               &pedge_vtx,
                                               PDM_OWNERSHIP_KEEP);

  /* Get faces */
  PDM_multipart_part_connectivity_get(mpart,
                                      0,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                      &pface_edge_idx,
                                      &pface_edge,
                                      PDM_OWNERSHIP_KEEP);

  /* Get groups */
  // _get_groups(mpart,
  //             0,
  //             i_part,
  //             &pn_edge_group,
  //             &pgroup_edge_idx,
  //             &pgroup_edge);
  PDM_g_num_t *edge_bound_ln_to_gn  = NULL;

  PDM_multipart_group_get(mpart, 0, i_part, PDM_MESH_ENTITY_EDGE,
                          &pn_edge_group,
                          &pgroup_edge_idx,
                          &pgroup_edge,
                          &edge_bound_ln_to_gn,
                          PDM_OWNERSHIP_KEEP);

  /* Transpose into edge_face */
  int   n_part          = 1;
  // int **tmp_pedge_face_idx = NULL;
  // int **tmp_pedge_face     = NULL;

  // PDM_part_connectivity_transpose(n_part,
  //                                 &pn_face,
  //                                 &pn_edge,
  //                                 &pface_edge_idx,
  //                                 &pface_edge,
  //                                 &tmp_pedge_face_idx,
  //                                 &tmp_pedge_face);

  // // int *pedge_face_idx = tmp_pedge_face_idx[i_part];
  // // int *pedge_face     = tmp_pedge_face[i_part];

  /* Transpose into edge_group */
  int **tmp_pedge_group_idx = NULL;
  int **tmp_pedge_group     = NULL;

  PDM_part_connectivity_transpose(n_part,
                                  &pn_edge_group,
                                  &pn_edge,
                                  &pgroup_edge_idx,
                                  &pgroup_edge,
                                  &tmp_pedge_group_idx,
                                  &tmp_pedge_group);

  int *pedge_group_idx = tmp_pedge_group_idx[i_part];
  int *pedge_group     = tmp_pedge_group    [i_part];
  free(tmp_pedge_group_idx);
  free(tmp_pedge_group    );

  /* Transpose into vtx_edge */
  int **tmp_pvtx_edge_idx = NULL;
  int **tmp_pvtx_edge     = NULL;

  int *pedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, pn_edge);

  PDM_part_connectivity_transpose(n_part,
                                  &pn_edge,
                                  &pn_vtx,
                                  &pedge_vtx_idx,
                                  &pedge_vtx,
                                  &tmp_pvtx_edge_idx,
                                  &tmp_pvtx_edge);

  int *pvtx_edge_idx = tmp_pvtx_edge_idx[i_part];
  int *pvtx_edge     = tmp_pvtx_edge    [i_part];
  free(tmp_pvtx_edge_idx);
  free(tmp_pvtx_edge    );

  /* Combine into vtx_group */

  int *pvtx_group_idx = NULL;
  int *pvtx_group     = NULL;

  PDM_combine_connectivity(pn_vtx,
                           pvtx_edge_idx,
                           pvtx_edge,
                           pedge_group_idx,
                           pedge_group,
                           &pvtx_group_idx,
                           &pvtx_group);

  /* Entities for Laplacian Smoothing */
  char filename[999];
  int                  pvtx_idx1             = 0;
  int                  pvtx_idx2             = 0;
  int                 *pstrid                = PDM_array_const_int (pn_vtx, 1);
  int                 *pnormalisation        = PDM_array_const_int (pn_vtx, 0);
  double              *pnew_vtx_coord        = malloc(3 * pn_vtx * sizeof(double));
  PDM_block_to_part_t *btp                   = NULL;
  PDM_part_to_block_t *ptb                   = NULL;
  int                 *dnormalisation        = NULL;
  int                 *dstrid_normalisation  = NULL;
  double              *dnew_vtx_coord        = NULL;
  int                 *dstrid_new_vtx_coord  = NULL;
  PDM_g_num_t         *distrib               = NULL;
  // PDM_g_num_t         *block_gnum            = NULL;
  int                  nelmt_proc            = 0;
  double              *dnormalisation_summed = NULL;
  double              *dnew_vtx_coord_summed = NULL;
  int                  idx                   = 0;
  int                  partial_sum           = 0;
  double               partial_sum_x         = 0;
  double               partial_sum_y         = 0;
  int                  stride                = 1;

  for (int i = 0; i < 3*pn_vtx; i++) {
      pnew_vtx_coord[i] = 0;
  }

  /* Laplacian smoothing */

  for (int i_step = 0; i_step <= n_steps; i_step++) {

    /* Export mesh to vtk format */
    if(1 == 0) {
      if (n_rank == 1) {
        sprintf(filename, "mesh_%2.2d.vtk", i_step);
      } else {
        sprintf(filename, "mesh_%2.2d_%2.2d.vtk", i_rank, i_step);
      }

      PDM_vtk_write_std_elements(filename,
                                 pn_vtx,
                                 pvtx_coord,
                                 vtx_ln_to_gn,
                                 PDM_MESH_NODAL_BAR2,
                                 pn_edge,
                                 pedge_vtx,
                                 NULL,
                                 0,
                                 NULL,
                                 NULL);
    }

    /* Local Laplacian Smoothing on partition */

    for (int i = 0; i < pn_edge; i++) {

      pvtx_idx1 = pedge_vtx[2*i]-1;
      pvtx_idx2 = pedge_vtx[2*i+1]-1;

      pnew_vtx_coord[3*pvtx_idx1]   += pvtx_coord[3*pvtx_idx2];
      pnew_vtx_coord[3*pvtx_idx1+1] += pvtx_coord[3*pvtx_idx2+1];
      pnew_vtx_coord[3*pvtx_idx2]   += pvtx_coord[3*pvtx_idx1];
      pnew_vtx_coord[3*pvtx_idx2+1] += pvtx_coord[3*pvtx_idx1+1];

      if (i_step == 0) {
        pnormalisation[pvtx_idx1]++;
        pnormalisation[pvtx_idx2]++;
      }

      // if (pvtx_group_idx[pvtx_idx1+1] - pvtx_group_idx[pvtx_idx1] == 0)

    } // end loop on partionned edges

    /* Merge coordinates using part_to_block */

    // Create part_to_block

    ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                   PDM_PART_TO_BLOCK_POST_MERGE,
                                   1.,
                                   &vtx_ln_to_gn,
                                   NULL,
                                   &pn_vtx,
                                   n_part,
                                   comm);

    // Exchange part_to_block

    if (i_step == 0) {

      PDM_part_to_block_exch(ptb,
                             sizeof(int),
                             PDM_STRIDE_VAR_INTERLACED,
                             1,
                             &pstrid,
                   (void **) &pnormalisation,
                             &dstrid_normalisation,
                   (void **) &dnormalisation);

    }

    int s_block_data = PDM_part_to_block_exch(ptb,
                           3 * sizeof(double),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           &pstrid,
                 (void **) &pnew_vtx_coord,
                           &dstrid_new_vtx_coord,
                 (void **) &dnew_vtx_coord);

    PDM_UNUSED (s_block_data);

    // Get block information

    distrib    = PDM_part_to_block_distrib_index_get(ptb);
    // block_gnum = PDM_part_to_block_block_gnum_get(ptb);
    nelmt_proc = PDM_part_to_block_n_elt_block_get(ptb);

    if (i_step ==0) {
      dnormalisation_summed = malloc(    nelmt_proc * sizeof(double));
      dnew_vtx_coord_summed = malloc(3 * nelmt_proc * sizeof(double));
    }


    // Sum merged normalisation and new_vtx_coord

    if (i_step == 0) {

      // for normalisation
      idx = 0;
      for (int i = 0; i < nelmt_proc; i++) {
        partial_sum = 0;
        for (int j = 0; j < dstrid_normalisation[i]; j++) {
          partial_sum += dnormalisation[idx++];
        }
        dnormalisation_summed[i] = 1. / (double) partial_sum;
        // log_trace("elmt #"PDM_FMT_G_NUM" : sum = %d, normalisation = %f\n", block_gnum[i], partial_sum, dnormalisation_summed[i]);
      }
    }

    // for new coordinates
    idx           = 0;
    for (int i = 0; i < nelmt_proc; i++) {
      partial_sum_x = 0;
      partial_sum_y = 0;
      // log_trace("point #"PDM_FMT_G_NUM" : stride = %d, normalisation = %f\n",
      //           block_gnum[i], dstrid_new_vtx_coord[i], dnormalisation_summed[i]);
      for (int j = 0; j < dstrid_new_vtx_coord[i]; j++) {
        // log_trace("  add %f %f * %d\n", dnew_vtx_coord[3*idx], dnew_vtx_coord[3*idx+1], dnormalisation[idx]);
        partial_sum_x += dnew_vtx_coord[3*idx];// * dnormalisation[idx];
        partial_sum_y += dnew_vtx_coord[3*idx+1];// * dnormalisation[idx];
        idx++;
      }
      dnew_vtx_coord_summed[3*i]   = partial_sum_x * dnormalisation_summed[i];
      dnew_vtx_coord_summed[3*i+1] = partial_sum_y * dnormalisation_summed[i];
      dnew_vtx_coord_summed[3*i+2] = 0.;
      // log_trace("  new_coord = %f %f\n", dnew_vtx_coord_summed[3*i], dnew_vtx_coord_summed[3*i+1]);
    }

    /* Back to result using block_to_part */

    // Create block_to_part

    btp = PDM_block_to_part_create(distrib,
           (const PDM_g_num_t  **) &vtx_ln_to_gn,
                                   &pn_vtx,
                                   n_part,
                                   comm);

    // Exchange in place block_to_part

    PDM_block_to_part_exch_in_place(btp,
                                    3 * sizeof(double),
                                    PDM_STRIDE_CST_INTERLACED,
                                    &stride,
                           (void *) dnew_vtx_coord_summed,
                                    &pstrid,
                          (void **) &pnew_vtx_coord);

    /* End loop Updates */
    for (int i = 0; i < pn_vtx; i++) {
      // Change pvtx_coord if has no group
      if (pvtx_group_idx[i+1] - pvtx_group_idx[i] == 0) {
        pvtx_coord[3*i]       = pnew_vtx_coord[3*i];
        pvtx_coord[3*i+1]     = pnew_vtx_coord[3*i+1];
        pnew_vtx_coord[3*i]   = 0;
        pnew_vtx_coord[3*i+1] = 0;
      }
    }


    PDM_part_to_block_free(ptb);
    PDM_block_to_part_free(btp);

    free(dstrid_new_vtx_coord);
    free(dnew_vtx_coord);

  } // end loop on Laplacian Smoothing steps


  free(pedge_vtx_idx    );
  free(pstrid        );
  free(pnormalisation);
  free(pnew_vtx_coord);
  free(dnormalisation_summed);
  free(dnew_vtx_coord_summed);
  free(dstrid_normalisation);
  free(dnormalisation);

  free(pvtx_edge_idx);
  free(pvtx_edge    );
  free(pvtx_group_idx);
  free(pvtx_group    );
  free(pedge_group_idx);
  free(pedge_group    );

  /* Free entities */
  // TO DO add free of mallocs
  PDM_multipart_free(mpart);

  PDM_MPI_Finalize();
  return 0;
}


