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
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_gnum.h"
#include "pdm_multipart.h"
#include "pdm_block_to_block.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_elmts_priv.h"
#include "pdm_predicate.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"


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
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

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
_read_args(int            argc,
           char         **argv,
           char         **filename,
           int           *dim,
           int           *visu)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }

    else if (strcmp(argv[i], "-d") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *dim = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static PDM_dmesh_nodal_t *
_read_gamma_mesh
(
 PDM_MPI_Comm   comm,
 const char    *filename,
 int            mesh_dim
 )
{
  int fix_orientation = 1;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dim = 0;


  PDM_g_num_t gn_vtx   = 0;
  PDM_g_num_t gn_edge  = 0;
  PDM_g_num_t gn_tria  = 0;
  PDM_g_num_t gn_tetra = 0;

  double      *gvtx_coord   = NULL;
  PDM_g_num_t *gedge_vtx    = NULL;
  int         *gedge_group  = NULL;
  PDM_g_num_t *gtria_vtx    = NULL;
  int         *gtria_group  = NULL;
  PDM_g_num_t *gtetra_vtx   = NULL;
  int         *gtetra_group = NULL;

  char line[999];

  if (i_rank == 0) {

    FILE *f = fopen(filename, "r");

    if (f == NULL) {
      PDM_error(__FILE__, __LINE__, 0, "Could not read file %s\n", filename);
    }

    while (1) {

      int stat = fscanf(f, "%s", line);

      if (stat == EOF) {
        // End of file
        break;
      }


      if (strstr(line, "Dimension") != NULL) {
        // Get dimension
        fscanf(f, "%d", &dim);
      }


      else if (strstr(line, "Vertices") != NULL) {
        // Get vertices
        long _gn_vtx;
        fscanf(f, "%ld", &_gn_vtx);
        gn_vtx = (PDM_g_num_t) _gn_vtx;

        gvtx_coord = malloc(sizeof(double) * gn_vtx * 3);
        for (PDM_g_num_t i = 0; i < gn_vtx; i++) {
          for (int j = 0; j < dim; j++) {
            fscanf(f, "%lf", &gvtx_coord[3*i + j]);
          }
          for (int j = dim; j < 3; j++) {
            gvtx_coord[3*i + j] = 0.;
          }

          int ref;
          fscanf(f, "%d", &ref);
        }
      }


      else if (strstr(line, "Edges") != NULL) {
        // Get edges
        long _gn_edge;
        fscanf(f, "%ld", &_gn_edge);
        gn_edge = (PDM_g_num_t) _gn_edge;

        gedge_vtx   = malloc(sizeof(PDM_g_num_t) * gn_edge * 2);
        gedge_group = malloc(sizeof(int        ) * gn_edge);
        for (PDM_g_num_t i = 0; i < gn_edge; i++) {
          for (int j = 0; j < 2; j++) {
            fscanf(f, PDM_FMT_G_NUM, &gedge_vtx[2*i + j]);
          }
          fscanf(f, "%d", &gedge_group[i]);
        }
      }


      else if (strstr(line, "Triangles") != NULL) {
        // Get triangles
        long _gn_tria;
        fscanf(f, "%ld", &_gn_tria);
        gn_tria = (PDM_g_num_t) _gn_tria;

        gtria_vtx   = malloc(sizeof(PDM_g_num_t) * gn_tria * 3);
        gtria_group = malloc(sizeof(int        ) * gn_tria);
        for (PDM_g_num_t i = 0; i < gn_tria; i++) {
          for (int j = 0; j < 3; j++) {
            fscanf(f, PDM_FMT_G_NUM, &gtria_vtx[3*i + j]);
          }
          fscanf(f, "%d", &gtria_group[i]);
        }
      }


      else if (strstr(line, "Tetrahedra") != NULL) {
        // Get tetrahedra
        long _gn_tetra;
        fscanf(f, "%ld", &_gn_tetra);
        gn_tetra = (PDM_g_num_t) _gn_tetra;

        gtetra_vtx   = malloc(sizeof(PDM_g_num_t) * gn_tetra * 4);
        gtetra_group = malloc(sizeof(int        ) * gn_tetra);
        for (PDM_g_num_t i = 0; i < gn_tetra; i++) {
          for (int j = 0; j < 4; j++) {
            fscanf(f, PDM_FMT_G_NUM, &gtetra_vtx[4*i + j]);
          }
          fscanf(f, "%d", &gtetra_group[i]);
        }
      }


    }
    fclose(f);


    if (fix_orientation) {

      int n_tria_flipped = 0;
      if (mesh_dim == 2) {
        for (PDM_g_num_t i = 0; i < gn_tria; i++) {
          PDM_g_num_t *tv = gtria_vtx + 3*i;

          double vol = PDM_predicate_orient2d(gvtx_coord + 3*(tv[0] - 1),
                                              gvtx_coord + 3*(tv[1] - 1),
                                              gvtx_coord + 3*(tv[2] - 1));
          if (vol < 0) {
            n_tria_flipped++;

            PDM_g_num_t tmp = tv[0];
            gtria_vtx[3*i  ] = tv[1];
            gtria_vtx[3*i+1] = tmp;
          }
        }
      }

      int n_tetra_flipped = 0;
      for (PDM_g_num_t i = 0; i < gn_tetra; i++) {
        PDM_g_num_t *tv = gtetra_vtx + 4*i;

        double vol = PDM_predicate_orient3d(gvtx_coord + 3*(tv[0] - 1),
                                            gvtx_coord + 3*(tv[1] - 1),
                                            gvtx_coord + 3*(tv[2] - 1),
                                            gvtx_coord + 3*(tv[3] - 1));
        if (vol < 0) {
          n_tetra_flipped++;

          PDM_g_num_t tmp = tv[0];
          gtetra_vtx[4*i  ] = tv[1];
          gtetra_vtx[4*i+1] = tmp;
        }
      }

      if (1) {
        printf("flipped %d triangles / "PDM_FMT_G_NUM"\n", n_tria_flipped, gn_tria);
        if (gn_tetra > 0) {
          printf("flipped %d tetrahedra / "PDM_FMT_G_NUM"\n", n_tetra_flipped, gn_tetra);
        }
      }
    }

    if (0) {
      log_trace("dim = %d\n", dim);
      log_trace("gn_vtx = "PDM_FMT_G_NUM"\n", gn_vtx);
      for (PDM_g_num_t i = 0; i < gn_vtx; i++) {
        log_trace("vtx "PDM_FMT_G_NUM" : %f %f %f\n",
                  i+1, gvtx_coord[3*i], gvtx_coord[3*i+1], gvtx_coord[3*i+2]);
      }


      log_trace("gn_edge = "PDM_FMT_G_NUM"\n", gn_edge);
      for (PDM_g_num_t i = 0; i < gn_edge; i++) {
        log_trace("edge "PDM_FMT_G_NUM" : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM", group %d\n",
                  i+1, gedge_vtx[2*i], gedge_vtx[2*i+1], gedge_group[i]);
      }


      log_trace("gn_tria = "PDM_FMT_G_NUM"\n", gn_tria);
      for (PDM_g_num_t i = 0; i < gn_tria; i++) {
        log_trace("tria "PDM_FMT_G_NUM" : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM", group %d\n",
                  i+1, gtria_vtx[3*i], gtria_vtx[3*i+1], gtria_vtx[3*i+2], gtria_group[i]);
      }

      log_trace("gn_tetra = "PDM_FMT_G_NUM"\n", gn_tetra);
      for (PDM_g_num_t i = 0; i < gn_tetra; i++) {
        log_trace("tetra "PDM_FMT_G_NUM" : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM", group %d\n",
                  i+1, gtetra_vtx[4*i], gtetra_vtx[4*i+1], gtetra_vtx[4*i+2], gtetra_vtx[4*i+3], gtetra_group[i]);
      }
    }
  }

  /* Distribute mesh */
  int dn_vtx   = (int) gn_vtx;
  int dn_edge  = (int) gn_edge;
  int dn_tria  = (int) gn_tria;
  int dn_tetra = (int) gn_tetra;
  PDM_g_num_t *init_distrib_vtx   = PDM_compute_entity_distribution(comm, dn_vtx);
  PDM_g_num_t *init_distrib_edge  = PDM_compute_entity_distribution(comm, dn_edge);
  PDM_g_num_t *init_distrib_tria  = PDM_compute_entity_distribution(comm, dn_tria);
  PDM_g_num_t *init_distrib_tetra = PDM_compute_entity_distribution(comm, dn_tetra);

  // PDM_MPI_Bcast(&dim,      1, PDM_MPI_INT,        0, comm);
  PDM_MPI_Bcast(&gn_vtx,   1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&gn_edge,  1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&gn_tria,  1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&gn_tetra, 1, PDM__PDM_MPI_G_NUM, 0, comm);

  /* Vertices*/
  PDM_g_num_t *distrib_vtx = PDM_compute_uniform_entity_distribution(comm, gn_vtx);
  PDM_block_to_block_t *btb_vtx = PDM_block_to_block_create(init_distrib_vtx,
                                                             distrib_vtx,
                                                             comm);
  double *dvtx_coord = NULL;
  PDM_block_to_block_exch(btb_vtx,
                          3*sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void *)  gvtx_coord,
                          NULL,
                (void **) &dvtx_coord);

  /* Edges */
  PDM_g_num_t *distrib_edge = PDM_compute_uniform_entity_distribution(comm, gn_edge);
  PDM_block_to_block_t *btb_edge = PDM_block_to_block_create(init_distrib_edge,
                                                             distrib_edge,
                                                             comm);
  PDM_g_num_t *dedge_vtx = NULL;
  PDM_block_to_block_exch(btb_edge,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          2,
                          NULL,
                (void *)  gedge_vtx,
                          NULL,
                (void **) &dedge_vtx);

  int *dedge_group = NULL;
  PDM_block_to_block_exch(btb_edge,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void *)  gedge_group,
                          NULL,
                (void **) &dedge_group);

  /* Triangles */
  PDM_g_num_t *distrib_tria = PDM_compute_uniform_entity_distribution(comm, gn_tria);

  PDM_block_to_block_t *btb_tria = PDM_block_to_block_create(init_distrib_tria,
                                                             distrib_tria,
                                                             comm);

  PDM_g_num_t *dtria_vtx = NULL;
  PDM_block_to_block_exch(btb_tria,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          3,
                          NULL,
                (void *)  gtria_vtx,
                          NULL,
                (void **) &dtria_vtx);

  int *dtria_group = NULL;
  PDM_block_to_block_exch(btb_tria,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void *)  gtria_group,
                          NULL,
                (void **) &dtria_group);

  /* Tetrahedra */
  PDM_g_num_t *distrib_tetra = PDM_compute_uniform_entity_distribution(comm, gn_tetra);

  PDM_block_to_block_t *btb_tetra = PDM_block_to_block_create(init_distrib_tetra,
                                                              distrib_tetra,
                                                              comm);

  PDM_g_num_t *dtetra_vtx = NULL;
  PDM_block_to_block_exch(btb_tetra,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          4,
                          NULL,
                (void *)  gtetra_vtx,
                          NULL,
                (void **) &dtetra_vtx);

  // int *dtetra_group = NULL;
  // PDM_block_to_block_exch(btb_tetra,
  //                         sizeof(int),
  //                         PDM_STRIDE_CST_INTERLACED,
  //                         1,
  //                         NULL,
  //               (void *)  gtetra_group,
  //                         NULL,
  //               (void **) &dtetra_group);


  if (gvtx_coord   != NULL) free(gvtx_coord);
  if (gedge_vtx    != NULL) free(gedge_vtx);
  if (gedge_group  != NULL) free(gedge_group);
  if (gtria_vtx    != NULL) free(gtria_vtx);
  if (gtria_group  != NULL) free(gtria_group);
  if (gtetra_vtx   != NULL) free(gtetra_vtx);
  if (gtetra_group != NULL) free(gtetra_group);

  PDM_block_to_block_free(btb_vtx);
  PDM_block_to_block_free(btb_edge);
  PDM_block_to_block_free(btb_tria);
  PDM_block_to_block_free(btb_tetra);

  free(init_distrib_vtx  );
  free(init_distrib_edge );
  free(init_distrib_tria );
  free(init_distrib_tetra);

  // int mesh_dimension = dim;
  // if (dim == 3 && gn_tetra == 0) {
  //   mesh_dimension = 2;
  // }

  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  mesh_dim,
                                                  gn_vtx,
                                                  gn_edge,
                                                  gn_tria,
                                                  gn_tetra);

  dn_vtx   = (int) (distrib_vtx[i_rank+1]   - distrib_vtx[i_rank]);
  dn_edge  = (int) (distrib_edge[i_rank+1]  - distrib_edge[i_rank]);
  dn_tria  = (int) (distrib_tria[i_rank+1]  - distrib_tria[i_rank]);
  dn_tetra = (int) (distrib_tetra[i_rank+1] - distrib_tetra[i_rank]);

  /* Vertices */
  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);


  /* Sections */
  dmn->ridge->n_g_elmts = gn_edge;
  int id_section_edge = PDM_DMesh_nodal_elmts_section_add(dmn->ridge,
                                                          PDM_MESH_NODAL_BAR2);
  PDM_DMesh_nodal_elmts_section_std_set(dmn->ridge,
                                        id_section_edge,
                                        dn_edge,
                                        dedge_vtx,
                                        PDM_OWNERSHIP_KEEP);



  dmn->surfacic->n_g_elmts = gn_tria;
  int id_section_tria = PDM_DMesh_nodal_elmts_section_add(dmn->surfacic,
                                                          PDM_MESH_NODAL_TRIA3);
  PDM_DMesh_nodal_elmts_section_std_set(dmn->surfacic,
                                        id_section_tria,
                                        dn_tria,
                                        dtria_vtx,
                                        PDM_OWNERSHIP_KEEP);

  if (mesh_dim == 3) {
    dmn->volumic->n_g_elmts = gn_tetra;
    int id_section_tetra = PDM_DMesh_nodal_elmts_section_add(dmn->volumic,
                                                             PDM_MESH_NODAL_TETRA4);
    PDM_DMesh_nodal_elmts_section_std_set(dmn->volumic,
                                          id_section_tetra,
                                          dn_tetra,
                                          dtetra_vtx,
                                          PDM_OWNERSHIP_KEEP);
  }


  /* Groups */
  int _n_group_edge = 0;
  for (int i = 0; i < dn_edge; i++) {
    _n_group_edge = PDM_MAX(_n_group_edge, dedge_group[i]);// refs are > 0 (?)
    dedge_group[i] -= 1;
  }

  int n_group_edge;
  PDM_MPI_Allreduce(&_n_group_edge, &n_group_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  // log_trace("n_group_edge = %d\n", n_group_edge);

  int *dedge_group_idx = PDM_array_new_idx_from_const_stride_int(1, dn_edge);
  // PDM_log_trace_connectivity_int(dedge_group_idx, dedge_group, dn_edge, "dedge_group : ");

  int         *dgroup_edge_idx = NULL;
  PDM_g_num_t *dgroup_edge     = NULL;
  PDM_dentity_group_transpose(n_group_edge,
                              dedge_group_idx,
                              dedge_group,
                              distrib_edge,
                              &dgroup_edge_idx,
                              &dgroup_edge,
                              dmn->comm);
  free(dedge_group_idx);
  free(dedge_group);

  PDM_DMesh_nodal_elmts_group_set(dmn->ridge,
                                  n_group_edge,
                                  dgroup_edge_idx,
                                  dgroup_edge,
                                  PDM_OWNERSHIP_KEEP);

  /* !!! At this point we assume all faces are tirangles (needs some changes if there are also quads) */
  int dn_face = dn_tria;// + dn_quad;
  PDM_g_num_t *distrib_face = distrib_tria;
  int *dface_group = dtria_group;
  int _n_group_face = 0;
  for (int i = 0; i < dn_tria; i++) {
    _n_group_face = PDM_MAX(_n_group_face, dface_group[i]);// refs are > 0 (?)
    dface_group[i] -= 1;
  }

  int n_group_face;
  PDM_MPI_Allreduce(&_n_group_face, &n_group_face, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  // log_trace("n_group_face = %d\n", n_group_face);

  int *dface_group_idx = PDM_array_new_idx_from_const_stride_int(1, dn_face);

  int         *dgroup_face_idx = NULL;
  PDM_g_num_t *dgroup_face     = NULL;
  PDM_dentity_group_transpose(n_group_face,
                              dface_group_idx,
                              dface_group,
                              distrib_face,
                              &dgroup_face_idx,
                              &dgroup_face,
                              dmn->comm);
  free(dface_group_idx);
  free(dtria_group);

  PDM_DMesh_nodal_elmts_group_set(dmn->surfacic,
                                  n_group_face,
                                  dgroup_face_idx,
                                  dgroup_face,
                                  PDM_OWNERSHIP_KEEP);

  free(distrib_vtx  );
  free(distrib_edge );
  free(distrib_tria );
  free(distrib_tetra);

  return dmn;
}






static void
_separate_groups
(
 PDM_dmesh_nodal_t   *dmn,
 PDM_g_num_t        **distrib_vtx,
 double             **dvtx_coord,
 int                 *n_ridge,
 PDM_g_num_t       ***ridge_distrib_edge,
 PDM_g_num_t       ***ridge_dedge_vtx,
 int                 *n_surface,
 PDM_g_num_t       ***surface_distrib_face,
 PDM_g_num_t       ***surface_dface_vtx,
 PDM_ownership_t      ownership
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  const PDM_g_num_t *_distrib_vtx = PDM_DMesh_nodal_distrib_vtx_get(dmn);

  int dn_vtx = (int) (_distrib_vtx[i_rank+1] - _distrib_vtx[i_rank]);
  double *_dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

  if (ownership == PDM_OWNERSHIP_USER) {
    *dvtx_coord = malloc(sizeof(double) * dn_vtx * 3);
    memcpy(*dvtx_coord, _dvtx_coord, sizeof(double) * dn_vtx * 3);

    *distrib_vtx = malloc(sizeof(PDM_g_num_t) * (n_rank+1));
    memcpy(*distrib_vtx, _distrib_vtx, sizeof(PDM_g_num_t) * (n_rank+1));
  } else {
    *dvtx_coord = _dvtx_coord;
    *distrib_vtx = (PDM_g_num_t *) _distrib_vtx;
  }


  /* Ridges */
  int         *dridge_edge_idx = NULL;
  PDM_g_num_t *dridge_edge     = NULL;
  PDM_DMesh_nodal_section_group_elmt_get(dmn,
                                         PDM_GEOMETRY_KIND_RIDGE,
                                         n_ridge,
                                         &dridge_edge_idx,
                                         &dridge_edge);


  int *ridge_dn_edge  = malloc(sizeof(int          ) * (*n_ridge));
  *ridge_distrib_edge = malloc(sizeof(PDM_g_num_t *) * (*n_ridge));
  for (int iridge = 0; iridge < *n_ridge; iridge++) {
    ridge_dn_edge[iridge] = dridge_edge_idx[iridge+1] - dridge_edge_idx[iridge];

    (*ridge_distrib_edge)[iridge] = PDM_compute_entity_distribution(dmn->comm,
                                                                    ridge_dn_edge[iridge]);
  }

  const PDM_g_num_t **pedge_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * (*n_ridge));
  for (int iridge = 0; iridge < *n_ridge; iridge++) {
    pedge_ln_to_gn[iridge] = dridge_edge + dridge_edge_idx[iridge];
  }

  {
    int *sections_id = PDM_DMesh_nodal_sections_id_get(dmn, PDM_GEOMETRY_KIND_RIDGE);
    int  n_section   = PDM_DMesh_nodal_n_section_get  (dmn, PDM_GEOMETRY_KIND_RIDGE);
    assert(n_section == 1);

    int id_section = sections_id[0];
    const PDM_g_num_t    *distrib_edge = PDM_DMesh_nodal_distrib_section_get(dmn, PDM_GEOMETRY_KIND_RIDGE, id_section);
    PDM_g_num_t          *dedge_vtx    = PDM_DMesh_nodal_section_std_get    (dmn, PDM_GEOMETRY_KIND_RIDGE, id_section);
    PDM_Mesh_nodal_elt_t  t_elt        = PDM_DMesh_nodal_section_type_get   (dmn, PDM_GEOMETRY_KIND_RIDGE, id_section);

    assert(t_elt == PDM_MESH_NODAL_BAR2 || t_elt == PDM_MESH_NODAL_BARHO);

    int order = 1;//?
    int edge_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

    PDM_block_to_part_t *btp_edge = PDM_block_to_part_create(distrib_edge,
                                                             pedge_ln_to_gn,
                                                             ridge_dn_edge,
                                                             *n_ridge,
                                                             dmn->comm);

    PDM_block_to_part_exch(btp_edge,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &edge_vtx_n,
                (void   *) dedge_vtx,
                           NULL,
                (void ***) ridge_dedge_vtx);
    PDM_block_to_part_free(btp_edge);
    free(pedge_ln_to_gn);
    free(ridge_dn_edge);
  }



  /* Surfaces */

  int         *dsurface_face_idx = NULL;
  PDM_g_num_t *dsurface_face     = NULL;
  PDM_DMesh_nodal_section_group_elmt_get(dmn,
                                         PDM_GEOMETRY_KIND_SURFACIC,
                                         n_surface,
                                         &dsurface_face_idx,
                                         &dsurface_face);


  int *surface_dn_face  = malloc(sizeof(int          ) * (*n_surface));
  *surface_distrib_face = malloc(sizeof(PDM_g_num_t *) * (*n_surface));
  for (int isurface = 0; isurface < *n_surface; isurface++) {
    surface_dn_face[isurface] = dsurface_face_idx[isurface+1] - dsurface_face_idx[isurface];

    (*surface_distrib_face)[isurface] = PDM_compute_entity_distribution(dmn->comm,
                                                                        surface_dn_face[isurface]);
  }

  const PDM_g_num_t **pface_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * (*n_surface));
  for (int isurface = 0; isurface < *n_surface; isurface++) {
    pface_ln_to_gn[isurface] = dsurface_face + dsurface_face_idx[isurface];
  }

  {
    int *sections_id = PDM_DMesh_nodal_sections_id_get(dmn, PDM_GEOMETRY_KIND_SURFACIC);
    int  n_section   = PDM_DMesh_nodal_n_section_get  (dmn, PDM_GEOMETRY_KIND_SURFACIC);
    assert(n_section == 1);

    int id_section = sections_id[0];
    const PDM_g_num_t    *distrib_face = PDM_DMesh_nodal_distrib_section_get(dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_g_num_t          *dface_vtx    = PDM_DMesh_nodal_section_std_get    (dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);
    PDM_Mesh_nodal_elt_t  t_elt        = PDM_DMesh_nodal_section_type_get   (dmn, PDM_GEOMETRY_KIND_SURFACIC, id_section);

    assert(t_elt == PDM_MESH_NODAL_TRIA3 || t_elt == PDM_MESH_NODAL_TRIAHO);

    int order = 1;//?
    int face_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

    PDM_block_to_part_t *btp_face = PDM_block_to_part_create(distrib_face,
                                                             pface_ln_to_gn,
                                                             surface_dn_face,
                                                             *n_surface,
                                                             dmn->comm);

    PDM_block_to_part_exch(btp_face,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST_INTERLACED,
                           &face_vtx_n,
                (void   *) dface_vtx,
                           NULL,
                (void ***) surface_dface_vtx);
    PDM_block_to_part_free(btp_face);
    free(pface_ln_to_gn);
    free(surface_dn_face);
  }

}




static void
_dump_dstd_elt
(
 const char           *prefix,
 PDM_MPI_Comm          comm,
 PDM_g_num_t          *distrib_vtx,
 double               *dvtx_coord,
 PDM_g_num_t          *distrib_elt,
 PDM_g_num_t          *delt_vtx,
 PDM_Mesh_nodal_elt_t  elt_type,
 const int             order
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  char filename[999];
  sprintf(filename, "%s_%2.2d.vtk", prefix, i_rank);

  int dn_elt = (int) (distrib_elt[i_rank+1] - distrib_elt[i_rank]);

  int elt_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, order);
  int *delt_vtx_idx = PDM_array_new_idx_from_const_stride_int(elt_vtx_n, dn_elt);

  PDM_g_num_t *pelt_ln_to_gn = malloc(sizeof(PDM_g_num_t) * dn_elt);
  for (int i = 0; i < dn_elt; i++) {
    pelt_ln_to_gn[i] = distrib_elt[i_rank] + i + 1;
  }

  PDM_g_num_t *pvtx_ln_to_gn;
  int         *pelt_vtx_idx;
  int         *pelt_vtx;
  int          pn_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_elt,
                                                           delt_vtx_idx,
                                                           delt_vtx,
                                                           dn_elt,
                                     (const PDM_g_num_t *) pelt_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pelt_vtx_idx,
                                                           &pelt_vtx);
  free(delt_vtx_idx);

  double **tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_vtx,
                                        dvtx_coord,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  double *pvtx_coord = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);

  assert(order == 1 && elt_type <= PDM_MESH_NODAL_HEXA8);

  PDM_vtk_write_std_elements(filename,
                             pn_vtx,
                             pvtx_coord,
                             pvtx_ln_to_gn,
                             elt_type,
                             dn_elt,
                             pelt_vtx,
                             pelt_ln_to_gn,
                             0,
                             NULL,
                             NULL);

  free(pelt_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pelt_vtx_idx);
  free(pelt_vtx);
  free(pvtx_coord);

}








int main(int argc, char *argv[])
{
  /*
   *  Init
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_predicate_exactinit();

  /*
   *  Read args
   */
  char *filename = NULL;
  int   mesh_dim = 3;
  int   visu     = 0;
  _read_args(argc,
             argv,
             &filename,
             &mesh_dim,
             &visu);

  if (filename == NULL) {
    filename = (char *) PDM_MESH_DIR"box.mesh";
  }


  /*
   *  Read initial mesh
   */
  PDM_dmesh_nodal_t *dmn = _read_gamma_mesh(comm,
                                            filename,
                                            mesh_dim);

  /*
   *  Separate groups
   */
  PDM_ownership_t ownership = PDM_OWNERSHIP_KEEP;

  double       *dvtx_coord           = NULL;
  int           n_ridge              = 0;
  PDM_g_num_t **ridge_distrib_edge   = NULL;
  PDM_g_num_t **ridge_dedge_vtx      = NULL;
  int           n_surface            = 0;
  PDM_g_num_t **surface_distrib_face = NULL;
  PDM_g_num_t **surface_dface_vtx    = NULL;
  PDM_g_num_t  *distrib_vtx          = NULL;
  _separate_groups(dmn,
                   &distrib_vtx,
                   &dvtx_coord,
                   &n_ridge,
                   &ridge_distrib_edge,
                   &ridge_dedge_vtx,
                   &n_surface,
                   &surface_distrib_face,
                   &surface_dface_vtx,
                   ownership);

  int order = 1;



  /*
   *  Visu
   */
  if (visu) {
    char prefix[999];

    for (int i = 0; i < n_ridge; i++) {

      sprintf(prefix, "ridge_%d", i);
      _dump_dstd_elt(prefix,
                     comm,
                     distrib_vtx,
                     dvtx_coord,
                     ridge_distrib_edge[i],
                     ridge_dedge_vtx[i],
                     PDM_MESH_NODAL_BAR2,
                     order);

    }


    for (int i = 0; i < n_surface; i++) {

      sprintf(prefix, "surface_%d", i);
      _dump_dstd_elt(prefix,
                     comm,
                     distrib_vtx,
                     dvtx_coord,
                     surface_distrib_face[i],
                     surface_dface_vtx[i],
                     PDM_MESH_NODAL_TRIA3,
                     order);

    }
  }


  /*
   *  Free memory
   */
  PDM_DMesh_nodal_free(dmn);

  for (int i = 0; i < n_ridge; i++) {
    free(ridge_dedge_vtx[i]);
    free(ridge_distrib_edge[i]);
  }
  free(ridge_dedge_vtx);
  free(ridge_distrib_edge);

  for (int i = 0; i < n_surface; i++) {
    free(surface_dface_vtx[i]);
    free(surface_distrib_face[i]);
  }
  free(surface_dface_vtx);
  free(surface_distrib_face);

  if (ownership == PDM_OWNERSHIP_USER) {
    free(dvtx_coord);
    free(distrib_vtx);
  }

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize ();

  return 0;
}
