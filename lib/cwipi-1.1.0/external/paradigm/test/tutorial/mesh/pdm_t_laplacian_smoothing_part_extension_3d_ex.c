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
#include "pdm_writer.h"


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
  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 3);

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

  for (int i = 0; i < 3*vtx_distrib[i_rank]; i++) {
    rand();
  }

  double noise = 0.49 * length / (double) (n_vtx_seg - 1);
  for (int i = 0; i < dn_vtx; i++) {
    double x = dvtx_coord[3*i];
    double y = dvtx_coord[3*i+1];
    double z = dvtx_coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      double dx = noise * _rand();
      if (x > 0. && x < length && y > 0 && y < length && z > 0 && z < length) {
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
_get_groups_and_bounds
(
 PDM_multipart_t     *multipart,
 const int            i_domain,
 const int            i_part,
       int           *n_face_group,
       int          **face_bound_idx,
       int          **face_bound,
       int          **face_part_bound_proc_idx,
       int          **face_part_bound_part_idx,
       int          **face_part_bound,
       int          **vtx_part_bound_proc_idx,
       int          **vtx_part_bound_part_idx,
       int          **vtx_part_bound,
       PDM_g_num_t  **face_bound_ln_to_gn
)
{
  PDM_multipart_group_get(multipart,
                          0,
                          i_part,
                          PDM_MESH_ENTITY_FACE,
                          n_face_group,
                          face_bound_idx,
                          face_bound,
                          face_bound_ln_to_gn,
                          PDM_OWNERSHIP_KEEP);

  PDM_multipart_part_graph_comm_get(multipart,
                                    i_domain,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    vtx_part_bound_proc_idx,
                                    vtx_part_bound_part_idx,
                                    vtx_part_bound,
                                    PDM_OWNERSHIP_KEEP);

  PDM_multipart_part_graph_comm_get(multipart,
                                    i_domain,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    face_part_bound_proc_idx,
                                    face_part_bound_part_idx,
                                    face_part_bound,
                                    PDM_OWNERSHIP_KEEP);

}



static void
_compute_face_vtx
(
 const int   n_face,
 int        *pface_edge_idx,
 int        *pface_edge,
 int        *pedge_vtx,
 int       **pface_vtx
 )
{
  int dbg = 0;
  *pface_vtx = (int *) malloc(sizeof(int) * pface_edge_idx[n_face]);

  for (int i = 0; i < n_face; i++) {

    if (dbg) {
      log_trace("\nFace %d\n", i);
      for (int idx_edge = pface_edge_idx[i]; idx_edge < pface_edge_idx[i+1]; idx_edge++) {
        int iedge = PDM_ABS(pface_edge[idx_edge]) - 1;
        log_trace("  edge %d: %d %d\n",
                  pface_edge[idx_edge],
                  pedge_vtx[2*iedge], pedge_vtx[2*iedge+1]);
      }
    }
    int *_pface_vtx = *pface_vtx + pface_edge_idx[i];

    int cur_vtx, next_vtx;
    int cur_edge = pface_edge[pface_edge_idx[i]];
    if (cur_edge < 0) {
      cur_edge = -cur_edge - 1;
      cur_vtx  = pedge_vtx[2*cur_edge+1];
      next_vtx = pedge_vtx[2*cur_edge  ];
    } else {
      cur_edge = cur_edge - 1;
      cur_vtx  = pedge_vtx[2*cur_edge  ];
      next_vtx = pedge_vtx[2*cur_edge+1];
    }

    for (int ivtx = 0; ivtx < pface_edge_idx[i+1] - pface_edge_idx[i]; ivtx++) {
      _pface_vtx[ivtx] = cur_vtx;

      for (int iedg = pface_edge_idx[i]; iedg < pface_edge_idx[i+1]; iedg++) {
        cur_edge = pface_edge[iedg];
        int vtx1, vtx2;
        if (cur_edge < 0) {
          cur_edge = -cur_edge - 1;
          vtx1 = pedge_vtx[2*cur_edge+1];
          vtx2 = pedge_vtx[2*cur_edge  ];
        } else {
          cur_edge = cur_edge - 1;
          vtx1 = pedge_vtx[2*cur_edge  ];
          vtx2 = pedge_vtx[2*cur_edge+1];
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
      for (int ivtx = 0; ivtx < pface_edge_idx[i+1] - pface_edge_idx[i]; ivtx++) {
        log_trace("%d ", _pface_vtx[ivtx]);
      }
      log_trace("\n");
    }

  }
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
   *  5 -> tetra
   *  6 -> pyramid
   *  7 -> prism
   *  8 -> hexa
   */

  PDM_g_num_t          n_vtx_seg = 5;  // Number of vtx on each side of the cube mesh
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TETRA4;  // Type of cells
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

  int          pn_cell             = 0;
  int          pn_face             = 0;
  int          pn_edge             = 0;
  int          pn_vtx              = 0;
  int         *pcell_face_idx      = NULL;
  int         *pcell_face          = NULL;
  int         *pface_edge_idx      = NULL;
  int         *pface_edge          = NULL;
  int         *pedge_vtx           = NULL;
  double      *pvtx_coord          = NULL;
  int          pn_face_group       = 0;
  int         *pgroup_face_idx     = NULL;
  int         *pgroup_face         = NULL;
  PDM_g_num_t *cell_ln_to_gn       = NULL;
  PDM_g_num_t *face_ln_to_gn       = NULL;
  PDM_g_num_t *edge_ln_to_gn       = NULL;
  PDM_g_num_t *vtx_ln_to_gn        = NULL;
  PDM_g_num_t *face_group_ln_to_gn = NULL;

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

  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_EDGE,
                                  &edge_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  /* Get faces */
  pn_face = PDM_multipart_part_connectivity_get(mpart,
                                               0,
                                               i_part,
                                               PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                               &pface_edge_idx,
                                               &pface_edge,
                                               PDM_OWNERSHIP_KEEP);

  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &face_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  /* Get cells */
  pn_cell = PDM_multipart_part_connectivity_get(mpart,
                                                0,
                                                i_part,
                                                PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                &pcell_face_idx,
                                                &pcell_face,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_CELL,
                                  &cell_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  /* Get groups and part bounds */
  int *face_part_bound_proc_idx = NULL;
  int *face_part_bound_part_idx = NULL;
  int *face_part_bound          = NULL;
  int *vtx_part_bound_proc_idx  = NULL;
  int *vtx_part_bound_part_idx  = NULL;
  int *vtx_part_bound           = NULL;
  _get_groups_and_bounds (mpart,
                          0,
                          i_part,
                          &pn_face_group,
                          &pgroup_face_idx,
                          &pgroup_face,
                          &face_part_bound_proc_idx,
                          &face_part_bound_part_idx,
                          &face_part_bound,
                          &vtx_part_bound_proc_idx,
                          &vtx_part_bound_part_idx,
                          &vtx_part_bound,
                          &face_group_ln_to_gn);

  printf("&vtx_part_bound_proc_idx: %ls\n", vtx_part_bound_proc_idx);

  /* Get face_group */
  int *pface_group_idx = NULL;
  int *pface_group     = NULL;

  PDM_connectivity_transpose(pn_face_group,
                             pn_face,
                             pgroup_face_idx,
                             pgroup_face,
                            &pface_group_idx,
                            &pface_group);

  /* Get face_vtx (TO DO get with wright order) */

  int *pedge_vtx_idx = PDM_array_new_idx_from_const_stride_int(2, pn_edge);
  // int *pface_vtx_idx = NULL;
  // int *pface_vtx     = NULL;

  // PDM_combine_connectivity(pn_face,
  //                          pface_edge_idx,
  //                          pface_edge,
  //                          pedge_vtx_idx,
  //                          pedge_vtx,
  //                          &pface_vtx_idx,
  //                          &pface_vtx);

  /* Get face_cell */

  int *pface_cell_idx = NULL;
  int *pface_cell     = NULL;

  PDM_connectivity_transpose(pn_cell,
                             pn_face,
                             pcell_face_idx,
                             pcell_face,
                             &pface_cell_idx,
                             &pface_cell);


  /* Get vtx_group */

  int *pedge_face_idx = NULL;
  int *pedge_face     = NULL;

  PDM_connectivity_transpose(pn_face,
                             pn_edge,
                             pface_edge_idx,
                             pface_edge,
                             &pedge_face_idx,
                             &pedge_face);

  int *pedge_group_idx = NULL;
  int *pedge_group     = NULL;

  PDM_combine_connectivity(pn_edge,
                           pedge_face_idx,
                           pedge_face,
                           pface_group_idx,
                           pface_group,
                           &pedge_group_idx,
                           &pedge_group);

  int *pvtx_edge_idx = NULL;
  int *pvtx_edge     = NULL;

  PDM_connectivity_transpose(pn_edge,
                             pn_vtx,
                             pedge_vtx_idx,
                             pedge_vtx,
                             &pvtx_edge_idx,
                             &pvtx_edge);

  int *pvtx_group_idx = NULL;
  int *pvtx_group     = NULL;

  PDM_combine_connectivity(pn_vtx,
                           pvtx_edge_idx,
                           pvtx_edge,
                           pedge_group_idx,
                           pedge_group,
                           &pvtx_group_idx,
                           &pvtx_group);



  // for (int i = 0; i < pn_vtx; i++) {
  //   log_trace("sommet d'indice %d à %d groupes qui sont:", i, (pvtx_group_idx[i+1]- pvtx_group_idx[i]));
  //   for (int j = pvtx_group_idx[i]; j < pvtx_group_idx[i+1]; j++) {
  //     log_trace(" %d ", pvtx_group[j]);
  //   }
  //   log_trace("\n");
  // }



  /* Build face_vtx */
  // int *pface_vtx_idx;
  int *pface_vtx;
  // get_ordered_face_vtx(&pface_vtx,
  //                      &pface_vtx_idx,
  //                      pface_edge_idx,
  //                      pface_edge,
  //                      pedge_vtx,
  //                      pn_face);

  _compute_face_vtx(pn_face,
                    pface_edge_idx,
                    pface_edge,
                    pedge_vtx,
                    &pface_vtx);

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
                              pn_cell,
                              pn_face,
                              0, // USELESS n_face_part_bound,
                              pn_face_group,
                              pn_edge,
                              pn_vtx,
                              pcell_face_idx,
                              pcell_face,
                              pface_cell,
                              pface_edge_idx,
                              pface_edge,
                              pface_edge_idx,
                              pface_vtx,
                              pedge_vtx,
                              pface_group_idx,
                              pface_group,
                              NULL, // face_join_idx
                              NULL, // face_join
                              face_part_bound_proc_idx,
                              face_part_bound_part_idx,
                              face_part_bound,
                              vtx_part_bound_proc_idx,
                              vtx_part_bound_part_idx,
                              vtx_part_bound,
                              cell_ln_to_gn,
                              face_ln_to_gn,
                              edge_ln_to_gn,
                              vtx_ln_to_gn,
                              face_group_ln_to_gn,
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

  // Get coordinates

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
  double *normalisation  = malloc(    pn_vtx * sizeof(double));
  double *pvtx_coord_new = malloc(3 * pn_vtx * sizeof(double));

  /* Set up of Ensight output */

  // PDM_writer_t *id_cs = PDM_writer_create("Ensight",
  //                                         PDM_WRITER_FMT_BIN,
  //                                         PDM_WRITER_TOPO_DEFORMABLE, // topologie constante mais coordonnées variables
  //                                         PDM_WRITER_OFF,
  //                                         "test_3d_ens",
  //                                         "lapalcian_smoothing",
  //                                         PDM_MPI_COMM_WORLD,
  //                                         PDM_IO_KIND_MPI_SIMPLE,
  //                                         1.,
  //                                         NULL);

  // int id_geom = PDM_writer_geom_create(id_cs,
  //                                      "test3d_geom",
  //                                      n_part);


  int *face_vtx_n  = (int *) malloc (sizeof(int) * pn_face);
  int *cell_face_n = (int *) malloc (sizeof(int) * pn_cell);

  for (int i = 0; i < pn_face; i++) {
    face_vtx_n[i] = pface_edge_idx[i+1] - pface_edge_idx[i];
  }

  for (int i = 0; i < pn_cell; i++) {
    cell_face_n[i] = pcell_face_idx[i+1] - pcell_face_idx[i];
  }

  // Step
  for (int i_step = 0; i_step <= n_steps; i_step++) {

    // Output in ensight format
    // PDM_writer_step_beg(id_cs, (double) i_step);

    // if (i_step == 0) {
    //   PDM_writer_geom_coord_set (id_cs,
    //                              id_geom,
    //                              0,
    //                              pn_vtx,
    //                              pvtx_coord,
    //                              vtx_ln_to_gn);
    //   PDM_writer_geom_cell3d_cellface_add (id_cs,
    //                                        id_geom,
    //                                        0,
    //                                        pn_cell,
    //                                        pn_face,
    //                                        pface_edge_idx,
    //                                        face_vtx_n,
    //                                        pface_vtx,
    //                                        pcell_face_idx,
    //                                        cell_face_n,
    //                                        pcell_face,
    //                                        cell_ln_to_gn);
    // }
    // PDM_writer_geom_write(id_cs,
    //                       id_geom);


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

    // Initialise pvtx_coord_new
    for (int i = 0; i < pn_vtx; i++) {
      pvtx_coord_new[3*i]   = 0;
      pvtx_coord_new[3*i+1] = 0;
      pvtx_coord_new[3*i+2] = 0;
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
      pvtx_coord_new[3*vtx1_idx+2] += pvtx_coord[3*vtx2_idx+2];
      pvtx_coord_new[3*vtx2_idx]   += pvtx_coord[3*vtx1_idx];
      pvtx_coord_new[3*vtx2_idx+1] += pvtx_coord[3*vtx1_idx+1];
      pvtx_coord_new[3*vtx2_idx+2] += pvtx_coord[3*vtx1_idx+2];

    } // end loop over own edges

    // Loop over extension edges

    for (int i = 0; i < pn_edge_extension; i++) {
      vtx1_idx = pedge_vtx_extension[2*i]-1;
      vtx2_idx = pedge_vtx_extension[2*i+1]-1;

      if ((vtx1_idx < pn_vtx) &&  (vtx2_idx >= pn_vtx)){
        if (i_step == 0) {
          normalisation[vtx1_idx]++;
        }

        pvtx_coord_new[3*vtx1_idx]   += pvtx_coord_extension[3*(vtx2_idx-pn_vtx)];
        pvtx_coord_new[3*vtx1_idx+1] += pvtx_coord_extension[3*(vtx2_idx-pn_vtx)+1];
        pvtx_coord_new[3*vtx1_idx+2] += pvtx_coord_extension[3*(vtx2_idx-pn_vtx)+2];
      }

      if ((vtx2_idx < pn_vtx) && (vtx1_idx >= pn_vtx)) {
        if (i_step == 0) {
          normalisation[vtx2_idx]++;
        }

        pvtx_coord_new[3*vtx2_idx]   += pvtx_coord_extension[3*(vtx1_idx-pn_vtx)];
        pvtx_coord_new[3*vtx2_idx+1] += pvtx_coord_extension[3*(vtx1_idx-pn_vtx)+1];
        pvtx_coord_new[3*vtx2_idx+2] += pvtx_coord_extension[3*(vtx1_idx-pn_vtx)+2];
      }

    } // end loop over extension edges

    // Add normalisation and update coordinates
    for (int i = 0; i < pn_vtx; i++) {
      if (i_step == 0) {
        normalisation[i] = 1 / normalisation[i];
      }

      if (pvtx_group_idx[i+1]- pvtx_group_idx[i] == 0) {
        pvtx_coord[3*i]   = pvtx_coord_new[3*i] * normalisation[i];
        pvtx_coord[3*i+1] = pvtx_coord_new[3*i+1] * normalisation[i];
        pvtx_coord[3*i+2] = pvtx_coord_new[3*i+2] * normalisation[i];
      }

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
      pvtx_coord_extension[3*i+2] = pvtx_coord_extension_new[i_part][3*i+2];
    } // end loop on extension coordinates

    free(pvtx_coord_extension_new[i_part]);
    free(pvtx_coord_extension_new);

    // PDM_writer_step_end(id_cs);
  } // end Laplacian Smoothing loop


  /* Free entities */
  PDM_multipart_free(mpart);
  PDM_part_extension_free(pe);
  PDM_part_to_part_free(ptp);

  free(face_vtx_n  );
  free(cell_face_n );
  free(extension_vtx_gnum_idx );

  free(pvtx_edge_idx);
  free(pvtx_edge);
  free(pedge_group_idx);
  free(pedge_group);
  free(pedge_face_idx);
  free(pedge_face);
  free(pedge_vtx_idx);
  free(pface_vtx);
  free(pface_cell_idx);
  free(pface_cell);
  free(pface_group_idx);
  free(pface_group);
  free(pvtx_group_idx);
  free(pvtx_group    );

  free(normalisation );
  free(pvtx_coord_new);

  PDM_MPI_Finalize();
  return 0;
}
