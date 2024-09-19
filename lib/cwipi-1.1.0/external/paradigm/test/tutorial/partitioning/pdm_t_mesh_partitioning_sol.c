#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_mesh_nodal.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_printf.h"
#include "pdm_part_extension.h"

#include "pdm_writer_priv.h"
#include "pdm_part_connectivity_transform.h"

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
     "  -vtk   Mesh vtk output.\n\n"
     "  -fe    Get FE mesh structure instead of FV.\n\n"
     "  -h     This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int    argc,
           char **argv,
           int   *vtk,
           int   *fe)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-vtk") == 0) {
      *vtk = 1;
    }

    else if (strcmp(argv[i], "-fe") == 0) {
      *fe = 1;
    }

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/**
 *
 * \brief  Main
 *
 */

int main
(
 int   argc,
 char *argv[]
)
{
  int vtk = 0;
  int fe  = 0;

  _read_args(argc,
             argv,
             &vtk,
             &fe);

  // Initialize MPI environment
  int          i_rank = -1;
  int          n_rank = -1;
  PDM_MPI_Comm comm   = PDM_MPI_COMM_WORLD;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Generate block-distributed parallelepided mesh
  int                  n_x      = 10;
  int                  n_y      = 10;
  int                  n_z      = 10;
  double               lengthx  = 1.;
  double               xmin     = 0.;
  double               ymin     = 0.;
  double               zmin     = 0.;
  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_TETRA4;
  int                  order    = 1; // call PDM_dcube_nodal_gen_ordering_set if order > 1
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_x,
                                                        n_y,
                                                        n_z,
                                                        lengthx,
                                                        xmin,
                                                        ymin,
                                                        zmin,
                                                        elt_type,
                                                        order,
                                                        PDM_OWNERSHIP_USER);

  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(dmn);

  // free
  PDM_dcube_nodal_gen_free(dcube);

  // Create partitioning object
  int              n_domain      = 1; // fixed
  int              n_part      = 1; // fixed
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH; // PDM_SPLIT_DUAL_WITH_HILBERT;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part,
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

  PDM_multipart_compute(mpart);

  int i_section = 0; // fixed
  int i_domain    = 0; // fixed
  int i_part    = 0; // fixed

  // Get mesh arrrays in FE structure
  if (fe) {
    PDM_part_mesh_nodal_t *pmn  = NULL;
    PDM_multipart_get_part_mesh_nodal(mpart,
                                      0,
                                      &pmn,
                                      PDM_OWNERSHIP_USER);

    int n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn,
                                                      i_section,
                                                      i_part);

    int         *connec              = NULL;
    PDM_g_num_t *numabs              = NULL;
    int         *parent_num          = NULL;
    PDM_g_num_t *parent_entity_g_num = NULL;
    PDM_part_mesh_nodal_section_std_get(pmn,
                                        i_section,
                                        i_part,
                                        &connec,
                                        &numabs,
                                        &parent_num,
                                        &parent_entity_g_num,
                                        PDM_OWNERSHIP_KEEP);

    int *elt_vtx_idx = PDM_array_new_idx_from_const_stride_int(4, n_elt);
    int *elt_vtx = malloc(sizeof(int) * elt_vtx_idx[n_elt]);
    memcpy(elt_vtx, connec, sizeof(int) * elt_vtx_idx[n_elt]);

    PDM_g_num_t *elt_ln_to_gn = PDM_part_mesh_nodal_g_num_get(pmn,
                                                              i_section,
                                                              i_part,
                                                              PDM_OWNERSHIP_USER);

    double *coords = NULL;
    int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                 i_domain,
                                                 i_part,
                                                 &coords,
                                                 PDM_OWNERSHIP_USER);

    PDM_g_num_t *vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmn,
                                                                  i_part);

    // Visualisation cell->vtx
    if (vtk) {
      char filename[999];
      sprintf(filename, "partitioning_FE_%3.3d.vtk", i_rank);
      PDM_vtk_write_std_elements(filename,
                                 n_vtx,
                                 coords,
                                 vtx_ln_to_gn,
                                 elt_type,
                                 n_elt,
                                 elt_vtx,
                                 elt_ln_to_gn,
                                 0,
                                 NULL,
                                 NULL);
    }

    // free
    free(elt_ln_to_gn);
    free(elt_vtx_idx);
    free(elt_vtx);
    free(coords);
    PDM_part_mesh_nodal_free(pmn);
  }

  // Get mesh arrrays in FV structure
  else {
    PDM_g_num_t *vtx_ln_to_gn = NULL;
    int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_domain,
                                                i_part,
                                                PDM_MESH_ENTITY_VTX,
                                                &vtx_ln_to_gn,
                                                PDM_OWNERSHIP_USER);

    double *coords = NULL;
    PDM_multipart_part_vtx_coord_get(mpart,
                                     i_domain,
                                     i_part,
                                     &coords,
                                     PDM_OWNERSHIP_USER);

    // PDM_g_num_t *edge_ln_to_gn = NULL;
    // int n_edge = PDM_multipart_part_ln_to_gn_get(mpart,
    //                                             i_domain,
    //                                             i_part,
    //                                             PDM_MESH_ENTITY_EDGE,
    //                                             &edge_ln_to_gn,
    //                                             PDM_OWNERSHIP_USER);

    // int *edge_vtx_idx = NULL;
    // int *edge_vtx     = NULL;
    // PDM_multipart_part_connectivity_get(mpart,
    //                                     i_domain,
    //                                     i_part,
    //                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
    //                                     &edge_vtx_idx,
    //                                     &edge_vtx,
    //                                     PDM_OWNERSHIP_USER);

    // if (edge_vtx_idx != NULL) free (edge_vtx_idx);

    PDM_g_num_t *face_ln_to_gn = NULL;
    int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_domain,
                                                i_part,
                                                PDM_MESH_ENTITY_FACE,
                                                &face_ln_to_gn,
                                                PDM_OWNERSHIP_USER);

    // int *face_edge_idx = NULL;
    // int *face_edge     = NULL;
    // PDM_multipart_part_connectivity_get(mpart,
    //                                     i_domain,
    //                                     i_part,
    //                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
    //                                     &face_edge_idx,
    //                                     &face_edge,
    //                                     PDM_OWNERSHIP_USER);

    int *face_vtx_idx = NULL;
    int *face_vtx     = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        i_domain,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &face_vtx_idx,
                                        &face_vtx,
                                        PDM_OWNERSHIP_USER);

    PDM_g_num_t *cell_ln_to_gn = NULL;
    int n_cell = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_domain,
                                                i_part,
                                                PDM_MESH_ENTITY_CELL,
                                                &cell_ln_to_gn,
                                                PDM_OWNERSHIP_USER);

    int *cell_face_idx = NULL;
    int *cell_face     = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        i_domain,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                        &cell_face_idx,
                                        &cell_face,
                                        PDM_OWNERSHIP_USER);

    // // Visualisation edge->vtx
    // if (vtk) {
    //   char filename[999];
    //   sprintf(filename, "partitioning_FV_%3.3d.vtk", i_rank);
    //   PDM_vtk_write_std_elements(filename,
    //                              n_vtx,
    //                              coords,
    //                              vtx_ln_to_gn,
    //                              PDM_MESH_NODAL_BAR2,
    //                              n_edge,
    //                              edge_vtx,
    //                              edge_ln_to_gn,
    //                              0,
    //                              NULL,
    //                              NULL);
    // }

    // Use PDM_compute_face_vtx_from_face_and_edge if you need face->vtx connectivity

    // BONUS

    // step 1 : create
    PDM_extend_type_t  extend_type = PDM_EXTEND_FROM_VTX;
    int                depth       = 1;
    PDM_part_extension_t *part_ext = PDM_part_extension_create(n_domain,
                                                               &n_part,
                                                               extend_type,
                                                               depth,
                                                               comm,
                                                               PDM_OWNERSHIP_KEEP);

    // step 2 : set
    int *vtx_part_bound_proc_idx = NULL;
    int *vtx_part_bound_part_idx = NULL;
    int *vtx_part_bound          = NULL;
    PDM_multipart_part_graph_comm_get(mpart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_part_bound_proc_idx,
                                      &vtx_part_bound_part_idx,
                                      &vtx_part_bound,
                                      PDM_OWNERSHIP_KEEP);

    // PDM_part_extension_set_part(part_ext,
    //                             i_domain,
    //                             i_part,
    //                             n_cell,
    //                             n_face,
    //                             0, // n_face_part_bound
    //                             0, // n_face_group
    //                             n_edge,
    //                             n_vtx,
    //                             cell_face_idx,
    //                             cell_face,
    //                             NULL, // face_cell
    //                             face_edge_idx,
    //                             face_edge,
    //                             NULL, // face_vtx_idx
    //                             NULL, // face_vtx
    //                             edge_vtx,
    //                             NULL, // face_group_idx
    //                             NULL, // face_group
    //                             NULL, // face_join_idx
    //                             NULL, // face_join
    //                             NULL, // face_part_bound_proc_idx
    //                             NULL, // face_part_bound_part_idx
    //                             NULL, // face_part_bound
    //                             vtx_part_bound_proc_idx,
    //                             vtx_part_bound_part_idx,
    //                             vtx_part_bound,
    //                             cell_ln_to_gn,
    //                             face_ln_to_gn,
    //                             edge_ln_to_gn,
    //                             vtx_ln_to_gn,
    //                             NULL, // face_group_ln_to_gn
    //                             coords);


    PDM_part_extension_connectivity_set(part_ext,
                                        i_domain,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                        cell_face_idx,
                                        cell_face);

    PDM_part_extension_connectivity_set(part_ext,
                                        i_domain,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        face_vtx_idx,
                                        face_vtx);

    PDM_part_extension_vtx_coord_set(part_ext,
                                     i_domain,
                                     i_part,
                                     coords);

    PDM_part_extension_ln_to_gn_set(part_ext,
                                    i_domain,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    n_cell,
                                    cell_ln_to_gn);

    PDM_part_extension_ln_to_gn_set(part_ext,
                                    i_domain,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    n_face,
                                    face_ln_to_gn);

    PDM_part_extension_ln_to_gn_set(part_ext,
                                    i_domain,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    n_vtx,
                                    vtx_ln_to_gn);

    PDM_part_extension_part_bound_graph_set(part_ext,
                                            i_domain,
                                            i_part,
                                            PDM_MESH_ENTITY_VTX,
                                            vtx_part_bound_proc_idx,
                                            vtx_part_bound_part_idx,
                                            vtx_part_bound);

    // step 3 : compute
    PDM_part_extension_compute(part_ext);

    // step 4 : get
    // Cell
    PDM_g_num_t *cell_ln_to_gn_ext = NULL;
    int n_cell_ext = PDM_part_extension_ln_to_gn_get (part_ext,
                                                      i_domain,
                                                      i_part,
                                                      PDM_MESH_ENTITY_CELL,
                                                      &cell_ln_to_gn_ext);
    PDM_UNUSED(n_cell_ext);

    int *cell_face_ext     = NULL;
    int *cell_face_ext_idx = NULL;
    PDM_part_extension_connectivity_get (part_ext,
                                         i_domain,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                         &cell_face_ext_idx,
                                         &cell_face_ext);

    // Face
    PDM_g_num_t *face_ln_to_gn_ext = NULL;
    int n_face_ext = PDM_part_extension_ln_to_gn_get (part_ext,
                                                      i_domain,
                                                      i_part,
                                                      PDM_MESH_ENTITY_FACE,
                                                      &face_ln_to_gn_ext);
    PDM_UNUSED(n_face_ext);

    // int *face_edge_ext     = NULL;
    // int *face_edge_ext_idx = NULL;
    // PDM_part_extension_connectivity_get (part_ext,
    //                                      i_domain,
    //                                      i_part,
    //                                      PDM_CONNECTIVITY_TYPE_FACE_EDGE,
    //                                      &face_edge_ext_idx,
    //                                      &face_edge_ext);
    int *face_vtx_ext     = NULL;
    int *face_vtx_ext_idx = NULL;
    PDM_part_extension_connectivity_get (part_ext,
                                         i_domain,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                         &face_vtx_ext_idx,
                                         &face_vtx_ext);

    // Edge
    // PDM_g_num_t *edge_ln_to_gn_ext = NULL;
    // int n_edge_ext = PDM_part_extension_ln_to_gn_get (part_ext,
    //                                                   i_domain,
    //                                                   i_part,
    //                                                   PDM_MESH_ENTITY_EDGE,
    //                                                   &edge_ln_to_gn_ext);

    // int *edge_vtx_ext     = NULL;
    // int *edge_vtx_ext_idx = NULL;
    // PDM_part_extension_connectivity_get (part_ext,
    //                                      i_domain,
    //                                      i_part,
    //                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
    //                                      &edge_vtx_ext_idx,
    //                                      &edge_vtx_ext);

    // Vertices
    PDM_g_num_t *vtx_ln_to_gn_ext = NULL;
    int n_vtx_ext = PDM_part_extension_ln_to_gn_get (part_ext,
                                                      i_domain,
                                                      i_part,
                                                      PDM_MESH_ENTITY_VTX,
                                                      &vtx_ln_to_gn_ext);
    PDM_UNUSED(n_vtx_ext);

    double *vtx_coord_ext = NULL;
    PDM_part_extension_vtx_coord_get(part_ext,
                                     i_domain,
                                     i_part,
                                     &vtx_coord_ext);

    // Fusion partition with partition extension

    int total_n_cell = n_cell + n_cell_ext;
    int total_n_face = n_face + n_face_ext;
    // int total_n_edge = n_edge + n_edge_ext;
    int total_n_vtx  = n_vtx  + n_vtx_ext;

    // Cell
    PDM_g_num_t *total_cell_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * total_n_cell);
    memcpy(total_cell_ln_to_gn,          cell_ln_to_gn,     sizeof(PDM_g_num_t) * n_cell);
    memcpy(total_cell_ln_to_gn + n_cell, cell_ln_to_gn_ext, sizeof(PDM_g_num_t) * n_cell_ext);

    int *total_cell_face_idx = malloc(sizeof(int) * (total_n_cell + 1));
    memcpy(total_cell_face_idx, cell_face_idx, sizeof(int) * (n_cell + 1));

    for (int i = 0; i <= n_cell_ext; i++) {
      total_cell_face_idx[n_cell + i] = cell_face_idx[n_cell] + cell_face_ext_idx[i];
    } // end loop on extension cells

    int *total_cell_face = malloc(sizeof(int) * total_cell_face_idx[total_n_cell]);
    memcpy(total_cell_face,                         cell_face,     sizeof(int) * cell_face_idx[n_cell]);
    memcpy(total_cell_face + cell_face_idx[n_cell], cell_face_ext, sizeof(int) * cell_face_ext_idx[n_cell_ext]);

    // Face
    PDM_g_num_t *total_face_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * total_n_face);
    memcpy(total_face_ln_to_gn,          face_ln_to_gn,     sizeof(PDM_g_num_t) * n_face);
    memcpy(total_face_ln_to_gn + n_face, face_ln_to_gn_ext, sizeof(PDM_g_num_t) * n_face_ext);

    // int *total_face_edge_idx = malloc(sizeof(int) * (total_n_face + 1));
    // memcpy(total_face_edge_idx, face_edge_idx, sizeof(int) * (n_face + 1));

    // for (int i = 0; i <= n_face_ext; i++) {
    //   total_face_edge_idx[n_face + i] = face_edge_idx[n_face] + face_edge_ext_idx[i];
    // } // end loop on extension faces

    // int *total_face_edge = malloc(sizeof(int) * total_face_edge_idx[total_n_face]);
    // memcpy(total_face_edge,                         face_edge,     sizeof(int) * face_edge_idx[n_face]);
    // memcpy(total_face_edge + face_edge_idx[n_face], face_edge_ext, sizeof(int) * face_edge_ext_idx[n_face_ext]);
    int *total_face_vtx_idx = malloc(sizeof(int) * (total_n_face + 1));
    memcpy(total_face_vtx_idx, face_vtx_idx, sizeof(int) * (n_face + 1));

    for (int i = 0; i <= n_face_ext; i++) {
      total_face_vtx_idx[n_face + i] = face_vtx_idx[n_face] + face_vtx_ext_idx[i];
    } // end loop on extension faces

    int *total_face_vtx = malloc(sizeof(int) * total_face_vtx_idx[total_n_face]);
    memcpy(total_face_vtx,                        face_vtx,     sizeof(int) * face_vtx_idx[n_face]);
    memcpy(total_face_vtx + face_vtx_idx[n_face], face_vtx_ext, sizeof(int) * face_vtx_ext_idx[n_face_ext]);

    // // Edge
    // PDM_g_num_t *total_edge_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * total_n_edge);
    // memcpy(total_edge_ln_to_gn,          edge_ln_to_gn,     sizeof(PDM_g_num_t) * n_edge);
    // memcpy(total_edge_ln_to_gn + n_edge, edge_ln_to_gn_ext, sizeof(PDM_g_num_t) * n_edge_ext);

    // int *total_edge_vtx = malloc(sizeof(int) * 2 * total_n_edge);
    // memcpy(total_edge_vtx,              edge_vtx,     sizeof(int) * 2 * n_edge);
    // memcpy(total_edge_vtx + 2 * n_edge, edge_vtx_ext, sizeof(int) * 2 * n_edge_ext);

    // Vtx
    PDM_g_num_t *total_vtx_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * total_n_vtx);
    memcpy(total_vtx_ln_to_gn,          vtx_ln_to_gn,     sizeof(PDM_g_num_t) * n_vtx);
    memcpy(total_vtx_ln_to_gn + n_vtx, vtx_ln_to_gn_ext, sizeof(PDM_g_num_t) * n_vtx_ext);

    double *total_coords = (double *) malloc(sizeof(double) * 3 * total_n_vtx);
    memcpy(total_coords,             coords,        sizeof(double) * 3 * n_vtx);
    memcpy(total_coords + 3 * n_vtx, vtx_coord_ext, sizeof(double) * 3 * n_vtx_ext);

    // // Create face->vtx
    // int *total_face_vtx = NULL;
    // PDM_compute_face_vtx_from_face_and_edge(total_n_face,
    //                                         total_face_edge_idx,
    //                                         total_face_edge,
    //                                         total_edge_vtx,
    //                                         &total_face_vtx);

    // int *total_face_vtx_idx = malloc(sizeof(int) * (total_n_face+1));
    // for (int i = 0; i < total_n_face + 1; i++) {
    //   total_face_vtx_idx[i] = 3 * i; // triangle
    // }

    double *total_cell_color = malloc(sizeof(double) * total_n_cell);
    for (int i = 0; i < n_cell; i++) {
      total_cell_color[i] = 2*i_rank;
    }
    for (int i = n_cell; i < total_n_cell; i++) {
      total_cell_color[i] = 2*i_rank+1;
    }

    const char  *field_name[]   = {"extension"};
    double     **field_value[1] = {&total_cell_color};

    writer_wrapper(comm,
                   "visu",
                   "pmesh",
                    1, // n_part
                    &total_n_vtx,
                    &total_coords,
                    &total_vtx_ln_to_gn,
                    &total_n_cell,
                    &total_face_vtx_idx,
                    &total_face_vtx,
                    &total_cell_ln_to_gn,
                    -1, // cell_t
                    &total_n_face,
                    &total_cell_face_idx,
                    &total_cell_face,
                    "Ensight",
                    1, // n_elt_field
                    field_name, // elt_field_name
                    field_value, // elt_field_values
                    0, // n_vtx_field
                    NULL, // vtx_field_name
                    NULL); // vtx_field_values

    // free fusion
    free(total_cell_color);
    free(total_cell_ln_to_gn);
    free(total_cell_face_idx);
    free(total_cell_face);
    free(total_face_ln_to_gn);
    // free(total_face_edge_idx);
    // free(total_face_edge);
    // free(total_edge_ln_to_gn);
    // free(total_edge_vtx);
    free(total_vtx_ln_to_gn);
    free(total_coords);
    free(total_face_vtx_idx);
    free(total_face_vtx);

    // step 5 : free
    PDM_part_extension_free(part_ext);

    // free
    free(vtx_ln_to_gn);
    free(coords);
    // free(edge_ln_to_gn);
    // free(edge_vtx_idx);
    // free(edge_vtx);
    free(face_ln_to_gn);
    // free(face_edge_idx);
    // free(face_edge);
    free(face_vtx_idx);
    free(face_vtx);
    free(cell_ln_to_gn);
    free(cell_face_idx);
    free(cell_face);
  }

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  // Finalize MPI environment
  PDM_MPI_Finalize();

  return 0;
}
