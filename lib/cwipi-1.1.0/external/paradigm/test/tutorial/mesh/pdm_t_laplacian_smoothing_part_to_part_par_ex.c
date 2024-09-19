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
#include "pdm_unique.h"
#include "pdm_part_to_part.h"

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
  #ifdef PDM_HAVE_PARMETIS
    PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_PARMETIS;
  #else
    PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_HILBERT;
  #endif

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

//   log_trace("n_face_group = %d\n", n_face_group);
//   log_trace("n_join_groups = %d\n", n_join_groups);


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

  int          pn_face         = 0;
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
  pn_face = PDM_multipart_part_connectivity_get(mpart,
                                               0,
                                               i_part,
                                               PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                               &pface_edge_idx,
                                               &pface_edge,
                                               PDM_OWNERSHIP_KEEP);
  PDM_UNUSED(pn_face);

  /* Get groups */
  // int          pn_edge_group        = 0;
  // int         *pedge_group          = NULL;
  // int         *pedge_group_idx      = NULL;
  PDM_g_num_t *edge_bound_ln_to_gn  = NULL;

  PDM_multipart_group_get(mpart, 0, i_part, PDM_MESH_ENTITY_EDGE,
                          &pn_edge_group,
                          &pgroup_edge_idx,
                          &pgroup_edge,
                          &edge_bound_ln_to_gn,
                          PDM_OWNERSHIP_KEEP);

  // _get_groups(mpart,
  //             0,
  //             i_part,
  //             &pn_edge_group,
  //             &pgroup_edge_idx,
  //             &pgroup_edge);

  // log_trace("pn_edge_group = %d\n", pn_edge_group);

  /* Create pvtx_vtx_gnum */

  // Transpose into vtx_edge
  int   n_part            = 1;
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
  int *pvtx_edge     = tmp_pvtx_edge[i_part];
  free(tmp_pvtx_edge_idx);
  free(tmp_pvtx_edge);

  // Combine into vtx_vtx

  int *pvtx_vtx_idx = NULL;
  int *pvtx_vtx     = NULL;

  PDM_combine_connectivity(pn_vtx,
                           pvtx_edge_idx,
                           pvtx_edge,
                           pedge_vtx_idx,
                           pedge_vtx,
                           &pvtx_vtx_idx,
                           &pvtx_vtx);
  free(pedge_vtx_idx);

  // Create vtx_vtx stride
  int *pstrid = malloc(pn_vtx * sizeof(int));
  for (int i = 0; i < pn_vtx; i++) {
    pstrid[i] = pvtx_vtx_idx[i+1] - pvtx_vtx_idx[i];
  }

  // Create pvtx_vtx_gnum
  int adapted_size = pvtx_vtx_idx[pn_vtx] - pn_vtx;
  int adapted_index = 0;
  int index_normal = 0;
  PDM_g_num_t *pvtx_vtx_gnum = malloc(adapted_size * sizeof(PDM_g_num_t));
  for (int i = 0; i < pn_vtx; i++) {
    for (int j = 0; j < pstrid[i]; j++) {
      if (vtx_ln_to_gn[pvtx_vtx[index_normal]-1] != vtx_ln_to_gn[i]) { // if not my gnum
        pvtx_vtx_gnum[adapted_index++] = vtx_ln_to_gn[pvtx_vtx[index_normal]-1];
      }
      index_normal++;
    }
    pstrid[i]--; // compensate removed self
  }

  // // Create pvtx_vtx_gnum
  // PDM_g_num_t *pvtx_vtx_gnum = malloc(pvtx_vtx_idx[pn_vtx] * sizeof(PDM_g_num_t));
  // for (int i = 0; i < pvtx_vtx_idx[pn_vtx]; i++) {
  //   pvtx_vtx_gnum[i] = vtx_ln_to_gn[pvtx_vtx[i]-1];
  // }

  if(0 == 1) {
    int indice1 = 0;
    for (int i = 0; i < pn_vtx; i++) {
      log_trace("point #"PDM_FMT_G_NUM " à pour voisins", vtx_ln_to_gn[i]);
      for (int j = 0; j < pstrid[i]; j++) {
        log_trace(" %d", pvtx_vtx_gnum[indice1++]);
      }
      log_trace("\n");
    }
  }

  // PDM_log_trace_array_long(pvtx_vtx_gnum, adapted_size, "pvtx_vtx_gnum : ");
  // PDM_log_trace_array_int(pstrid, pn_vtx, "pstrid : ");

  /* Export mesh to vtk format */
  char filename[999];
  if(0 == 1) {
    sprintf(filename, "mesh_%2.2d.vtk", i_rank);

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
  /* part_to_block */

  int *dstrid_vtx_vtx_gnum  = NULL;
  PDM_g_num_t *dvtx_vtx_gnum        = NULL;

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &vtx_ln_to_gn,
                                                      NULL,
                                                      &pn_vtx,
                                                      n_part,
                                                      comm);

  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &pstrid,
               (void **) &pvtx_vtx_gnum,
                         &dstrid_vtx_vtx_gnum,
               (void **) &dvtx_vtx_gnum);

  int nelmt_proc = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *distrib = PDM_part_to_block_distrib_index_get(ptb);

  /* Remove duplicate */

  int *dstrid_vtx_vtx_gnum_sorted = malloc(nelmt_proc *sizeof(int));
  int idx_comp  = 0; // Compressed index use to fill the buffer
  int idx_block = 0; // Index in the block to post-treat

  int val = 0;
  int *tmp_idx = PDM_array_new_idx_from_sizes_int(dstrid_vtx_vtx_gnum, nelmt_proc);
  // PDM_log_trace_connectivity_long(tmp_idx, dvtx_vtx_gnum, nelmt_proc, "before dvtx_vtx_gnum : ");
  free(tmp_idx);

  for (int i = 0; i < nelmt_proc; i++) {

    val = PDM_inplace_unique_long(dvtx_vtx_gnum, NULL, idx_comp, idx_comp + dstrid_vtx_vtx_gnum[i]-1);
    dstrid_vtx_vtx_gnum_sorted[i] = val;
    idx_block += dstrid_vtx_vtx_gnum[i];
    idx_comp  += val;

    // (dstrid_vtx_vtx_gnum[i] != dstrid_vtx_vtx_gnum_sorted[i]) && (i+1 != nelmt_proc)
    if (i+1 != nelmt_proc){ // not last block
      for (int j = 0; j < dstrid_vtx_vtx_gnum[i+1]; j++) {
        dvtx_vtx_gnum[idx_comp + j] = dvtx_vtx_gnum[idx_block + j]; // shift next block
      } // end loop on neighbours to be shifted
    } // end if shifted
  } // end loop on vertices
  free(dstrid_vtx_vtx_gnum);

  // PDM_log_trace_array_int(dstrid_vtx_vtx_gnum_sorted, nelmt_proc, "dstrid_vtx_vtx_gnum_sorted: ");
  // PDM_log_trace_array_long(dvtx_vtx_gnum, idx_comp, "dvtx_vtx_gnum : ");
  tmp_idx = PDM_array_new_idx_from_sizes_int(dstrid_vtx_vtx_gnum_sorted, nelmt_proc);
  // PDM_log_trace_connectivity_long(tmp_idx, dvtx_vtx_gnum, nelmt_proc, "after dvtx_vtx_gnum : ");
  free(tmp_idx);
  free(pvtx_vtx_gnum);

  /* block_to_part */
  int         **pstrid_new       = NULL;
  PDM_g_num_t **pvtx_vtx_gnum_new = NULL;

  PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib,
                              (const PDM_g_num_t  **) &vtx_ln_to_gn,
                                                      &pn_vtx,
                                                      n_part,
                                                      comm);

  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         dstrid_vtx_vtx_gnum_sorted,
                (void *) dvtx_vtx_gnum,
                         &pstrid_new,
              (void ***) &pvtx_vtx_gnum_new);
  free(dstrid_vtx_vtx_gnum_sorted);
  free(dvtx_vtx_gnum);

  int size = PDM_block_to_part_n_elt_get(btp, i_part);
  int size_neighbours = 0;
  for (int i = 0; i < size; i++) {
    size_neighbours += pstrid_new[i_part][i];
  }
  PDM_UNUSED(size_neighbours);

  // int indice = 0;
  // for (int i = 0; i < size; i++) {
  //   log_trace("point #"PDM_FMT_G_NUM " à pour voisins", vtx_ln_to_gn[i]);
  //   for (int j = 0; j < pstrid_new[i_part][i]; j++) {
  //     log_trace(" %d", pvtx_vtx_gnum_new[i_part][indice++]);
  //   }
  // log_trace("\n");
  // }

  /* Laplacian Smoothing using part_to_part */

  // Fix boundaries

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
  free(tmp_pedge_group);

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
  free(pvtx_edge_idx);
  free(pvtx_edge    );
  // part_to_part

  // Create

  int *pvtx_vtx_gnum_new_idx = PDM_array_new_idx_from_sizes_int(pstrid_new[i_part], size);

  PDM_part_to_part_t *ptp = PDM_part_to_part_create((const PDM_g_num_t **) &vtx_ln_to_gn,
                                                                           &pn_vtx,
                                                                           1,
                                                    (const PDM_g_num_t **) &vtx_ln_to_gn,
                                                                           &pn_vtx,
                                                                           1,
                                                    (const int **)         &pvtx_vtx_gnum_new_idx,
                                                    (const PDM_g_num_t **) &pvtx_vtx_gnum_new[i_part],
                                                                           comm);

  // Needed entities
  int      request;
  double **pvtx_coord_neighbours = NULL;
  double  *tmp_coord             = malloc(3 * sizeof(double));

  // Step
  for (int i_step = 0; i_step <= n_steps; i_step++) {
    // Output in vtk format

    if(0 == 1) {
      sprintf(filename, "mesh_%2.2d_%2.2d.vtk", i_rank, i_step);
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

    // Get coordinates

    PDM_part_to_part_reverse_iexch(ptp,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   1,
                                   3 * sizeof(double),
                                   NULL,
                   (const void **) &pvtx_coord,
                                   NULL,
                        (void ***) &pvtx_coord_neighbours,
                                   &request);

    PDM_part_to_part_reverse_iexch_wait(ptp, request);

    // Compute new coordinates

    for (int i = 0; i < pn_vtx; i++) {
      tmp_coord[0] = 0;
      tmp_coord[1] = 0;
      for (int j = pvtx_vtx_gnum_new_idx[i]; j < pvtx_vtx_gnum_new_idx[i+1]; j++) {
        tmp_coord[0] += pvtx_coord_neighbours[i_part][3*j];
        tmp_coord[1] += pvtx_coord_neighbours[i_part][3*j+1];
      } // end loop on neighbour vetices
      // Fix boundaries
      if (pvtx_group_idx[i+1] - pvtx_group_idx[i] == 0) {
        pvtx_coord[3*i] = tmp_coord[0] / pstrid_new[i_part][i]; // faire normalisation tab
        pvtx_coord[3*i+1] = tmp_coord[1] / pstrid_new[i_part][i];
      } // end if no group
    } // end loop on vertices

    free(pvtx_coord_neighbours[i_part]);
    free(pvtx_coord_neighbours);
  } // end loop Laplace Smoothing stepping

  /* Free entities */

  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);
  PDM_part_to_part_free (ptp);
  PDM_multipart_free(mpart);
  free(tmp_coord);
  free(pstrid_new[i_part]);
  free(pstrid_new);
  free(pstrid);
  free(pvtx_vtx_gnum_new[i_part]);
  free(pvtx_vtx_gnum_new);
  free(pvtx_vtx_idx);
  free(pvtx_vtx);
  free(pedge_group_idx);
  free(pedge_group    );
  free(pvtx_group_idx );
  free(pvtx_group     );
  free(pvtx_vtx_gnum_new_idx);

  PDM_MPI_Finalize();
  return 0;
}


