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
#include "pdm_part_connectivity_transform.h"
#include "pdm_unique.h"
#include "pdm_part_to_part.h"
#include "pdm_part_extension.h"


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




static void
_get_groups
(
 PDM_multipart_t     *multipart,
 const int            i_domain,
 const int            i_part,
       int           *n_face_group,
       int          **face_bound_idx,
       int          **face_bound
 )
{
  PDM_g_num_t *face_bound_ln_to_gn = NULL;
  PDM_multipart_group_get(multipart,
                          0,
                          i_part,
                          PDM_MESH_ENTITY_FACE,
                          n_face_group,
                          face_bound_idx,
                          face_bound,
                          &face_bound_ln_to_gn,
                          PDM_OWNERSHIP_KEEP);

}




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
  PDM_Mesh_nodal_elt_t elt_type  = 2;  // Type of cells
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

  /* Get groups */
  _get_groups(mpart,
              0,
              i_part,
              &pn_edge_group,
              &pgroup_edge_idx,
              &pgroup_edge);

  /* TO DO get face_vtx */

  /* TO DO get vtx_group */

  /* part_extension */

  // Create

  int n_part = 1;

  PDM_part_extension_t *pe =  PDM_part_extension_create(1,
                                                        &n_part,
                                                        PDM_EXTEND_FROM_VTX,
                                                        1, // extension depth
                                                        comm,
                                                        PDM_OWNERSHIP_KEEP);

  // Set

  PDM_part_extension_set_part(pe,
                              0, // i_domain
                              i_part,
                              0, // n_cell
                              pn_face,
                              0, // n_face_part_bound
                              0, // n_face_group
                              pn_edge,
                              pn_vtx,
                              NULL, // cell_face_idx
                              NULL, // cell_face
                              NULL, // face_cell
                              pface_edge_idx,
                              pface_edge,
                              NULL, // face_vtx_idx
                              NULL, // face_vtx
                              pedge_vtx,
                              NULL, // face_bound_idx
                              NULL, // face_bound
                              NULL, // face_join_idx
                              NULL, // face_join
                              NULL, // face_part_bound_proc_idx
                              NULL, // face_part_bound_part_idx
                              NULL, // face_part_bound
                              NULL, // vtx_part_bound_proc_idx
                              NULL, // vtx_part_bound_part_idx
                              NULL, // vtx_part_bound
                              NULL, // cell_ln_to_gn
                              NULL, // face_ln_to_gn
                              NULL, // edge_ln_to_gn
                              vtx_ln_to_gn,
                              NULL, // face_group_ln_to_gn
                              pvtx_coord);

  // Compute

  PDM_part_extension_compute(pe);

  // Get extension connectivity

  int  pn_edge_extension;
  int *pedge_vtx_extension = NULL;
  int *pedge_vtx_extension_idx = NULL;

  pn_edge_extension = PDM_part_extension_connectivity_get(pe,
                                                          0, // i_domain
                                                          i_part,
                                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                          &pedge_vtx_extension_idx,
                                                          &pedge_vtx_extension);

  // Get coordinates (TO DO PDM_part_extension_group_get for extension ou group)

  double *pvtx_coord_extension = NULL;

  PDM_part_extension_vtx_coord_get(pe,
                               0, // i_domain
                               i_part,
                               &pvtx_coord_extension);

  // Get extension_vtx_gnum

  int          pn_vtx_extension;
  PDM_g_num_t *extension_vtx_gnum = NULL;

  pn_vtx_extension = PDM_part_extension_ln_to_gn_get(pe,
                                                     0, // i_domain
                                                     i_part,
                                                     PDM_MESH_ENTITY_VTX,
                                                     &extension_vtx_gnum);

  /* part_to_part */

  // Create

  int  request;
  int *extension_vtx_gnum_idx = PDM_array_new_idx_from_const_stride_int(1, pn_vtx_extension);

  PDM_part_to_part_t *ptp = PDM_part_to_part_create((const PDM_g_num_t **) &extension_vtx_gnum, // extension_vtx_gnum
                                                                           &pn_vtx_extension,
                                                                           1, // n_part1
                                                    (const PDM_g_num_t **) &vtx_ln_to_gn,
                                                                           &pn_vtx,
                                                                           1, // n_part2
                                                    (const int **)         &extension_vtx_gnum_idx, // stride 1
                                                    (const PDM_g_num_t **) &extension_vtx_gnum, // extension_vtx_gnum
                                                                           comm);

  /* Laplacian Smoothing */

  char    filename[999];
  int     vtx1_idx;
  int     vtx2_idx;
  double *normalisation  = malloc(pn_vtx * sizeof(double));
  double *pvtx_coord_new = malloc(3 * pn_vtx * sizeof(double));

  // Step
  for (int i_step = 0; i_step <= n_steps; i_step++) {

    // Output in vtk format (TO DO face_vtx)

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

    // Initialise pvtx_coord_new
    for (int i = 0; i < pn_vtx; i++) {
      pvtx_coord_new[i] = 0;
    }

    // Loop over own edges

    for (int i = 0; i < pn_edge; i++) {
      vtx1_idx = pedge_vtx[2*i]-1;
      vtx2_idx = pedge_vtx[2*i+1]-1;

      if (i_step == 0) {
        normalisation[vtx1_idx]++;
        normalisation[vtx2_idx]++;
      }

      pvtx_coord_new[3*vtx1_idx]   += pvtx_coord[3*vtx2_idx];
      pvtx_coord_new[3*vtx1_idx+1] += pvtx_coord[3*vtx2_idx+1];
      pvtx_coord_new[3*vtx2_idx]   += pvtx_coord[3*vtx1_idx];
      pvtx_coord_new[3*vtx2_idx+1] += pvtx_coord[3*vtx1_idx+1];

    } // end loop over own edges

    // Loop over extension edges

    for (int i = 0; i < pn_edge_extension; i++) {
      vtx1_idx = pedge_vtx[2*i]-1;
      vtx2_idx = pedge_vtx[2*i+1]-1;

      if (vtx1_idx <= pn_vtx) {
        if (i_step == 0) {
          normalisation[vtx1_idx]++;
        }

        pvtx_coord_new[3*vtx1_idx]   += pvtx_coord_extension[3*(vtx2_idx-pn_vtx)];
        pvtx_coord_new[3*vtx1_idx+1] += pvtx_coord_extension[3*(vtx2_idx-pn_vtx)+1];
      }

      if (vtx2_idx <= pn_vtx) {
        if (i_step == 0) {
          normalisation[vtx2_idx]++;
        }

        pvtx_coord_new[3*vtx2_idx]   += pvtx_coord_extension[3*(vtx1_idx-pn_vtx)];
        pvtx_coord_new[3*vtx2_idx+1] += pvtx_coord_extension[3*(vtx1_idx-pn_vtx)+1];
      }

    } // end loop over extension edges

    // Add normalisation and update coordinates (TO DO vtx_group to fix boundary)
    for (int i = 0; i < pn_edge; i++) {
      if (i_step == 0) {
        normalisation[i] = 1 / normalisation[i];
      }

      pvtx_coord[3*i]   = pvtx_coord_new[3*i] * normalisation[i];
      pvtx_coord[3*i+1] = pvtx_coord_new[3*i+1] * normalisation[i];

    } // end loop on coordinates

    // Update coordinates of extension vertices

    double **pvtx_coord_extension_new = NULL;

    PDM_part_to_part_reverse_iexch(ptp,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   1,
                                   3 * sizeof(double),
                                   NULL,
                   (const void **) &pvtx_coord,
                                   NULL,
                        (void ***) &pvtx_coord_extension_new,
                                   &request);

    PDM_part_to_part_reverse_iexch_wait(ptp, request);

    for (int i = 0; i < pn_vtx_extension; i++) {
      pvtx_coord_extension[3*i]   = pvtx_coord_extension_new[i_part][3*i];
      pvtx_coord_extension[3*i+1] = pvtx_coord_extension_new[i_part][3*i+1];
    } // end loop on extension coordinates

    free(pvtx_coord_extension_new);

  } // end Laplacian Smoothing loop

  /* Free entities */

  PDM_multipart_free(mpart);
  PDM_part_extension_free(pe);
  PDM_part_to_part_free(ptp);

  PDM_MPI_Finalize();
  return 0;
}
