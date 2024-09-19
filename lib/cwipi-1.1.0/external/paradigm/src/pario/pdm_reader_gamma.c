/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <sys/types.h>
#include <limits.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_block_to_block.h"
#include "pdm_part_to_block.h"
#include "pdm_distrib.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_elmts_priv.h"
#include "pdm_predicate.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

#include "pdm_reader_gamma.h"

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
 * Private function prototypes
 *============================================================================*/

static void _shift_groups
(
 const int  n_elt,
       int *elt_group
 )
{
  int min_group = INT_MAX;

  for (int i = 0; i < n_elt; i++) {
    min_group = PDM_MIN(min_group, elt_group[i]);
  }

  for (int i = 0; i < n_elt; i++) {
    elt_group[i] += 1 - min_group;
  }
}


/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a dmesh nodal from a file in ASCII GAMMA mesh format
 *
 * \param[in]  comm                MPI communicator
 * \param[in]  filename            Filename
 * \param[in]  fix_orientation_2d  Ensure positive area for 2d faces
 * \param[in]  fix_orientation_3d  Ensure positive volume for 3d cells
 *
 * \return Pointer to PDM_dmesh_nodal object
 *
 */

PDM_dmesh_nodal_t *
PDM_reader_gamma_dmesh_nodal
(
 PDM_MPI_Comm   comm,
 const char    *filename,
 int            fix_orientation_2d,
 int            fix_orientation_3d
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int dim = 0;


  PDM_g_num_t gn_vtx   = 0;
  PDM_g_num_t gn_edge  = 0;
  PDM_g_num_t gn_tria  = 0;
  PDM_g_num_t gn_tetra = 0;
  PDM_g_num_t gn_pyra  = 0;
  PDM_g_num_t gn_prism = 0;
  PDM_g_num_t gn_hexa  = 0;

  double      *gvtx_coord   = NULL;
  int         *gvtx_tag     = NULL;
  PDM_g_num_t *gedge_vtx    = NULL;
  int         *gedge_group  = NULL;
  PDM_g_num_t *gtria_vtx    = NULL;
  int         *gtria_group  = NULL;
  PDM_g_num_t *gtetra_vtx   = NULL;
  int         *gtetra_group = NULL;
  PDM_g_num_t *gpyra_vtx    = NULL;
  int         *gpyra_group  = NULL;
  PDM_g_num_t *gprism_vtx   = NULL;
  int         *gprism_group = NULL;
  PDM_g_num_t *ghexa_vtx    = NULL;
  int         *ghexa_group  = NULL;

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
        gvtx_tag   = malloc(sizeof(int   ) * gn_vtx);
        for (PDM_g_num_t i = 0; i < gn_vtx; i++) {
          for (int j = 0; j < dim; j++) {
            fscanf(f, "%lf", &gvtx_coord[3*i + j]);
          }
          for (int j = dim; j < 3; j++) {
            gvtx_coord[3*i + j] = 0.;
          }
          fscanf(f, "%d", &gvtx_tag[i]);
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


      else if (strstr(line, "Pyramids") != NULL) {
        // Get pyramids
        long _gn_pyra;
        fscanf(f, "%ld", &_gn_pyra);
        gn_pyra = (PDM_g_num_t) _gn_pyra;

        gpyra_vtx   = malloc(sizeof(PDM_g_num_t) * gn_pyra * 5);
        gpyra_group = malloc(sizeof(int        ) * gn_pyra);
        for (PDM_g_num_t i = 0; i < gn_pyra; i++) {
          for (int j = 0; j < 5; j++) {
            fscanf(f, PDM_FMT_G_NUM, &gpyra_vtx[5*i + j]);
          }
          fscanf(f, "%d", &gpyra_group[i]);
        }
      }


      else if (strstr(line, "Prisms") != NULL) {
        // Get prisms
        long _gn_prism;
        fscanf(f, "%ld", &_gn_prism);
        gn_prism = (PDM_g_num_t) _gn_prism;

        gprism_vtx   = malloc(sizeof(PDM_g_num_t) * gn_prism * 6);
        gprism_group = malloc(sizeof(int        ) * gn_prism);
        for (PDM_g_num_t i = 0; i < gn_prism; i++) {
          for (int j = 0; j < 6; j++) {
            fscanf(f, PDM_FMT_G_NUM, &gprism_vtx[6*i + j]);
          }
          fscanf(f, "%d", &gprism_group[i]);
        }
      }


      else if (strstr(line, "Hexahedra") != NULL) {
        // Get hexahedra
        long _gn_hexa;
        fscanf(f, "%ld", &_gn_hexa);
        gn_hexa = (PDM_g_num_t) _gn_hexa;

        ghexa_vtx   = malloc(sizeof(PDM_g_num_t) * gn_hexa * 8);
        ghexa_group = malloc(sizeof(int        ) * gn_hexa);
        for (PDM_g_num_t i = 0; i < gn_hexa; i++) {
          for (int j = 0; j < 8; j++) {
            fscanf(f, PDM_FMT_G_NUM, &ghexa_vtx[8*i + j]);
          }
          fscanf(f, "%d", &ghexa_group[i]);
        }
      }


    }
    fclose(f);


    if (fix_orientation_2d) {

      int n_tria_flipped = 0;
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

      if (0) {
        printf("flipped %d triangles / "PDM_FMT_G_NUM"\n", n_tria_flipped, gn_tria);
      }
    }

    if (gn_tetra > 0 && fix_orientation_3d) {
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

      if (0) {
        printf("flipped %d tetrahedra / "PDM_FMT_G_NUM"\n", n_tetra_flipped, gn_tetra);
      }
    }

    if (1) {
      /* Shift groups */
      // _shift_groups((int) gn_vtx,   gvtx_group); // TODO: when corners
      _shift_groups((int) gn_edge,  gedge_group);
      _shift_groups((int) gn_tria,  gtria_group);
      _shift_groups((int) gn_tetra, gtetra_group);
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
  int dn_pyra  = (int) gn_pyra;
  int dn_prism = (int) gn_prism;
  int dn_hexa  = (int) gn_hexa;
  PDM_g_num_t *init_distrib_vtx   = PDM_compute_entity_distribution(comm, dn_vtx);
  PDM_g_num_t *init_distrib_edge  = PDM_compute_entity_distribution(comm, dn_edge);
  PDM_g_num_t *init_distrib_tria  = PDM_compute_entity_distribution(comm, dn_tria);
  PDM_g_num_t *init_distrib_tetra = PDM_compute_entity_distribution(comm, dn_tetra);
  PDM_g_num_t *init_distrib_pyra  = PDM_compute_entity_distribution(comm, dn_pyra);
  PDM_g_num_t *init_distrib_prism = PDM_compute_entity_distribution(comm, dn_prism);
  PDM_g_num_t *init_distrib_hexa  = PDM_compute_entity_distribution(comm, dn_hexa);

  int mesh_dimension = dim;
  if (dim == 3 && gn_tetra == 0 && gn_pyra == 0 && gn_prism == 0 && gn_hexa == 0) {
    mesh_dimension = 2;
  }

  PDM_MPI_Bcast(&gn_vtx,   1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&gn_edge,  1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&gn_tria,  1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&gn_tetra, 1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&gn_pyra,  1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&gn_prism, 1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&gn_hexa,  1, PDM__PDM_MPI_G_NUM, 0, comm);
  PDM_MPI_Bcast(&mesh_dimension, 1, PDM_MPI_INT, 0, comm);

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

  int *dvtx_tag = NULL;
  PDM_block_to_block_exch(btb_vtx,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void *)  gvtx_tag,
                          NULL,
                (void **) &dvtx_tag);


  /* Edges */
  PDM_g_num_t *distrib_edge = PDM_compute_uniform_entity_distribution(comm, gn_edge);
  PDM_block_to_block_t *btb_edge = PDM_block_to_block_create(init_distrib_edge,
                                                             distrib_edge,
                                                             comm);
  PDM_g_num_t *dedge_vtx = NULL;
  if (gn_edge > 0) {
    PDM_block_to_block_exch(btb_edge,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            2,
                            NULL,
                  (void  *) gedge_vtx,
                            NULL,
                  (void **) &dedge_vtx);
  }

  int *dedge_group = NULL;
  PDM_block_to_block_exch(btb_edge,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void  *) gedge_group,
                          NULL,
                (void **) &dedge_group);

  /* Triangles */
  PDM_g_num_t *distrib_tria = PDM_compute_uniform_entity_distribution(comm, gn_tria);

  PDM_block_to_block_t *btb_tria = PDM_block_to_block_create(init_distrib_tria,
                                                             distrib_tria,
                                                             comm);

  PDM_g_num_t *dtria_vtx = NULL;
  if (gn_tria > 0) {
    PDM_block_to_block_exch(btb_tria,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            3,
                            NULL,
                  (void *)  gtria_vtx,
                            NULL,
                  (void **) &dtria_vtx);
  }

  int *dtria_group = NULL;
  PDM_block_to_block_exch(btb_tria,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void  *) gtria_group,
                          NULL,
                (void **) &dtria_group);

  /* Tetrahedra */
  PDM_g_num_t *distrib_tetra = PDM_compute_uniform_entity_distribution(comm, gn_tetra);

  PDM_block_to_block_t *btb_tetra = PDM_block_to_block_create(init_distrib_tetra,
                                                              distrib_tetra,
                                                              comm);

  PDM_g_num_t *dtetra_vtx = NULL;
  if (gn_tetra > 0) {
    // Can overflow in buffer send/recv - Two much stride
    // PDM_block_to_block_exch(btb_tetra,
    //                         sizeof(PDM_g_num_t),
    //                         PDM_STRIDE_CST_INTERLACED,
    //                         4,
    //                         NULL,
    //               (void  *) gtetra_vtx,
    //                         NULL,
    //               (void **) &dtetra_vtx);

    // Economic exchange
    PDM_MPI_Datatype mpi_tetra_type;
    PDM_MPI_Type_create_contiguous(4, PDM__PDM_MPI_G_NUM, &mpi_tetra_type);
    PDM_MPI_Type_commit(&mpi_tetra_type);

    PDM_block_to_block_exch_with_mpi_type(btb_tetra,
                                          PDM_STRIDE_CST_INTERLACED,
                                          mpi_tetra_type,
                                          NULL,
                                (void  *) gtetra_vtx,
                                          NULL,
                                (void **) &dtetra_vtx);

    PDM_MPI_Type_free(&mpi_tetra_type);
  }

  int *dtetra_group = NULL;
  PDM_block_to_block_exch(btb_tetra,
                          sizeof(int),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
                (void *)  gtetra_group,
                          NULL,
                (void **) &dtetra_group);


  /* Pyramids */
  PDM_g_num_t *distrib_pyra = PDM_compute_uniform_entity_distribution(comm, gn_pyra);

  PDM_block_to_block_t *btb_pyra = PDM_block_to_block_create(init_distrib_pyra,
                                                              distrib_pyra,
                                                              comm);

  PDM_g_num_t *dpyra_vtx = NULL;
  if (gn_pyra > 0) {
    PDM_block_to_block_exch(btb_pyra,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            5,
                            NULL,
                  (void  *) gpyra_vtx,
                            NULL,
                  (void **) &dpyra_vtx);
  }

  // int *dpyra_group = NULL;
  // PDM_block_to_block_exch(btb_pyra,
  //                         sizeof(int),
  //                         PDM_STRIDE_CST_INTERLACED,
  //                         1,
  //                         NULL,
  //               (void *)  gpyra_group,
  //                         NULL,
  //               (void **) &dpyra_group);


  /* Prisms */
  PDM_g_num_t *distrib_prism = PDM_compute_uniform_entity_distribution(comm, gn_prism);

  PDM_block_to_block_t *btb_prism = PDM_block_to_block_create(init_distrib_prism,
                                                              distrib_prism,
                                                              comm);

  PDM_g_num_t *dprism_vtx = NULL;
  if (gn_prism > 0) {
    PDM_block_to_block_exch(btb_prism,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            6,
                            NULL,
                  (void  *) gprism_vtx,
                            NULL,
                  (void **) &dprism_vtx);
  }

  // int *dprism_group = NULL;
  // PDM_block_to_block_exch(btb_prism,
  //                         sizeof(int),
  //                         PDM_STRIDE_CST_INTERLACED,
  //                         1,
  //                         NULL,
  //               (void *)  gprism_group,
  //                         NULL,
  //               (void **) &dprism_group);


  /* Hexahedra */
  PDM_g_num_t *distrib_hexa = PDM_compute_uniform_entity_distribution(comm, gn_hexa);

  PDM_block_to_block_t *btb_hexa = PDM_block_to_block_create(init_distrib_hexa,
                                                              distrib_hexa,
                                                              comm);

  PDM_g_num_t *dhexa_vtx = NULL;
  if (gn_hexa > 0) {
    PDM_block_to_block_exch(btb_hexa,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST_INTERLACED,
                            8,
                            NULL,
                  (void  *) ghexa_vtx,
                            NULL,
                  (void **) &dhexa_vtx);
  }

  // int *dhexa_group = NULL;
  // PDM_block_to_block_exch(btb_hexa,
  //                         sizeof(int),
  //                         PDM_STRIDE_CST_INTERLACED,
  //                         1,
  //                         NULL,
  //               (void *)  ghexa_group,
  //                         NULL,
  //               (void **) &dhexa_group);


  if (gvtx_coord   != NULL) free(gvtx_coord);
  if (gvtx_tag     != NULL) free(gvtx_tag  );
  if (gedge_vtx    != NULL) free(gedge_vtx);
  if (gedge_group  != NULL) free(gedge_group);
  if (gtria_vtx    != NULL) free(gtria_vtx);
  if (gtria_group  != NULL) free(gtria_group);
  if (gtetra_vtx   != NULL) free(gtetra_vtx);
  if (gtetra_group != NULL) free(gtetra_group);
  if (gpyra_vtx    != NULL) free(gpyra_vtx);
  if (gpyra_group  != NULL) free(gpyra_group);
  if (gprism_vtx   != NULL) free(gprism_vtx);
  if (gprism_group != NULL) free(gprism_group);
  if (ghexa_vtx    != NULL) free(ghexa_vtx);
  if (ghexa_group  != NULL) free(ghexa_group);

  PDM_block_to_block_free(btb_vtx);
  PDM_block_to_block_free(btb_edge);
  PDM_block_to_block_free(btb_tria);
  PDM_block_to_block_free(btb_tetra);
  PDM_block_to_block_free(btb_pyra);
  PDM_block_to_block_free(btb_prism);
  PDM_block_to_block_free(btb_hexa);

  free(init_distrib_vtx  );
  free(init_distrib_edge );
  free(init_distrib_tria );
  free(init_distrib_tetra);
  free(init_distrib_pyra );
  free(init_distrib_prism);
  free(init_distrib_hexa );

  if (i_rank == 0) {
    if (0) {
      printf("mesh_dimension = %d\n", mesh_dimension);
      printf("gn_vtx   = "PDM_FMT_G_NUM"\n", gn_vtx);
      printf("gn_edge  = "PDM_FMT_G_NUM"\n", gn_edge);
      printf("gn_tria  = "PDM_FMT_G_NUM"\n", gn_tria);
      printf("gn_tetra = "PDM_FMT_G_NUM"\n", gn_tetra);
      printf("gn_pyra  = "PDM_FMT_G_NUM"\n", gn_pyra);
      printf("gn_prism = "PDM_FMT_G_NUM"\n", gn_prism);
      printf("gn_hexa  = "PDM_FMT_G_NUM"\n", gn_hexa);
    }
  }

  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  mesh_dimension,
                                                  gn_vtx,
                                                  gn_tetra + gn_pyra + gn_prism + gn_hexa,
                                                  gn_tria,
                                                  gn_edge);

  dn_vtx   = (int) (distrib_vtx  [i_rank+1] - distrib_vtx  [i_rank]);
  dn_edge  = (int) (distrib_edge [i_rank+1] - distrib_edge [i_rank]);
  dn_tria  = (int) (distrib_tria [i_rank+1] - distrib_tria [i_rank]);
  dn_tetra = (int) (distrib_tetra[i_rank+1] - distrib_tetra[i_rank]);
  dn_pyra  = (int) (distrib_pyra [i_rank+1] - distrib_pyra [i_rank]);
  dn_prism = (int) (distrib_prism[i_rank+1] - distrib_prism[i_rank]);
  dn_hexa  = (int) (distrib_hexa [i_rank+1] - distrib_hexa [i_rank]);

  /* Vertices */
  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  /* Sections */
  dmn->ridge->n_g_elmts = gn_edge;
  if (gn_edge > 0) {
    int id_section_edge = PDM_DMesh_nodal_elmts_section_add(dmn->ridge,
                                                            PDM_MESH_NODAL_BAR2);
    PDM_DMesh_nodal_elmts_section_std_set(dmn->ridge,
                                          id_section_edge,
                                          dn_edge,
                                          dedge_vtx,
                                          PDM_OWNERSHIP_KEEP);
  } else {
    free(dedge_vtx);
  }


  dmn->surfacic->n_g_elmts = gn_tria;
  if (gn_tria > 0) {
    int id_section_tria = PDM_DMesh_nodal_elmts_section_add(dmn->surfacic,
                                                            PDM_MESH_NODAL_TRIA3);
    PDM_DMesh_nodal_elmts_section_std_set(dmn->surfacic,
                                          id_section_tria,
                                          dn_tria,
                                          dtria_vtx,
                                          PDM_OWNERSHIP_KEEP);
  } else {
    free(dtria_vtx);
  }

  if (mesh_dimension == 3) {
    dmn->volumic->n_g_elmts = gn_tetra + gn_pyra + gn_prism + gn_hexa;

    if (gn_tetra > 0) {
      int id_section_tetra = PDM_DMesh_nodal_elmts_section_add(dmn->volumic,
                                                               PDM_MESH_NODAL_TETRA4);
      PDM_DMesh_nodal_elmts_section_std_set(dmn->volumic,
                                            id_section_tetra,
                                            dn_tetra,
                                            dtetra_vtx,
                                            PDM_OWNERSHIP_KEEP);
    } else {
      free(dtetra_vtx);
    }

    if (gn_pyra > 0) {
      int id_section_pyra = PDM_DMesh_nodal_elmts_section_add(dmn->volumic,
                                                              PDM_MESH_NODAL_PYRAMID5);
      PDM_DMesh_nodal_elmts_section_std_set(dmn->volumic,
                                            id_section_pyra,
                                            dn_pyra,
                                            dpyra_vtx,
                                            PDM_OWNERSHIP_KEEP);
    } else {
      free(dpyra_vtx);
    }

    if (gn_prism > 0) {
      int id_section_prism = PDM_DMesh_nodal_elmts_section_add(dmn->volumic,
                                                               PDM_MESH_NODAL_PRISM6);
      PDM_DMesh_nodal_elmts_section_std_set(dmn->volumic,
                                            id_section_prism,
                                            dn_prism,
                                            dprism_vtx,
                                            PDM_OWNERSHIP_KEEP);
    } else {
      free(dprism_vtx);
    }

    if (gn_hexa > 0) {
      int id_section_hexa = PDM_DMesh_nodal_elmts_section_add(dmn->volumic,
                                                              PDM_MESH_NODAL_HEXA8);
      PDM_DMesh_nodal_elmts_section_std_set(dmn->volumic,
                                            id_section_hexa,
                                            dn_hexa,
                                            dhexa_vtx,
                                            PDM_OWNERSHIP_KEEP);
    } else {
      free(dhexa_vtx);
    }
  }


  /* Set vtx tags */
  PDM_DMesh_nodal_vtx_tag_set(dmn, dvtx_tag);

  /* Edge groups */
  int _n_group_edge = 0;
  for (int i = 0; i < dn_edge; i++) {
    _n_group_edge = PDM_MAX(_n_group_edge, dedge_group[i]);// refs are > 0 (?)
    dedge_group[i] -= 1;
  }
  // PDM_log_trace_array_int(dedge_group, dn_edge, "dedge_group : ");

  int n_group_edge;
  PDM_MPI_Allreduce(&_n_group_edge, &n_group_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  int *dedge_group_idx = PDM_array_new_idx_from_const_stride_int(1, dn_edge);

  int         *dgroup_edge_idx = NULL;
  PDM_g_num_t *dgroup_edge     = NULL;
  if(n_group_edge > 0) {
    PDM_dentity_group_transpose(n_group_edge,
                              dedge_group_idx,
                              dedge_group,
                              distrib_edge,
                              &dgroup_edge_idx,
                              &dgroup_edge,
                              dmn->comm);
  }
  free(dedge_group_idx);
  free(dedge_group);

  // PDM_log_trace_array_int(dgroup_edge_idx, n_group_edge+1, "dgroup_edge_idx : ");

  PDM_DMesh_nodal_elmts_group_set(dmn->ridge,
                                  n_group_edge,
                                  dgroup_edge_idx,
                                  dgroup_edge,
                                  PDM_OWNERSHIP_KEEP);

  /* Tria groups */
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

  // PDM_log_trace_array_int(dgroup_face_idx, n_group_face+1, "dgroup_face_idx : ");

  PDM_DMesh_nodal_elmts_group_set(dmn->surfacic,
                                  n_group_face,
                                  dgroup_face_idx,
                                  dgroup_face,
                                  PDM_OWNERSHIP_KEEP);


  /* Tetra groups */
  int _n_group_tetra = 0;
  for (int i = 0; i < dn_tetra; i++) {
    _n_group_tetra = PDM_MAX(_n_group_tetra, dtetra_group[i]);// refs are > 0 (?)
    dtetra_group[i] -= 1;
  }

  int n_group_tetra;
  PDM_MPI_Allreduce(&_n_group_tetra, &n_group_tetra, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  int *dtetra_group_idx = PDM_array_new_idx_from_const_stride_int(1, dn_tetra);

  int         *dgroup_tetra_idx = NULL;
  PDM_g_num_t *dgroup_tetra     = NULL;
  PDM_dentity_group_transpose(n_group_tetra,
                              dtetra_group_idx,
                              dtetra_group,
                              distrib_tetra,
                              &dgroup_tetra_idx,
                              &dgroup_tetra,
                              dmn->comm);
  free(dtetra_group_idx);
  free(dtetra_group);

  // PDM_log_trace_array_int(dgroup_face_idx, n_group_face+1, "dgroup_face_idx : ");

  PDM_DMesh_nodal_elmts_group_set(dmn->volumic,
                                  n_group_tetra,
                                  dgroup_tetra_idx,
                                  dgroup_tetra,
                                  PDM_OWNERSHIP_KEEP);

  free(distrib_vtx  );
  free(distrib_edge );
  free(distrib_tria );
  free(distrib_tetra);
  free(distrib_pyra );
  free(distrib_prism);
  free(distrib_hexa );

  return dmn;
}



void
PDM_write_meshb
(
  const char   *filename,
  const int     n_vtx,
  const int     n_tetra,
  const int     n_tri,
  const int     n_edge,
  const double *vtx_coords,
  const int    *vtx_tags,
  const int    *tetra_vtx,
  const int    *tetra_tag,
  const int    *tria_vtx,
  const int    *tria_tag,
  const int    *edge_vtx,
  const int    *edge_tag
)
{
  // Write file
  FILE *f = fopen(filename, "w");

  fprintf(f, "MeshVersionFormatted 2\n");
  fprintf(f, "# rank %d\n\n", 0);
  fprintf(f, "Dimension\n3\n\n");

  fprintf(f, "Vertices\n%d\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    fprintf(f, "%20.16lf %20.16lf %20.16lf %i\n",
            vtx_coords[3*i  ],
            vtx_coords[3*i+1],
            vtx_coords[3*i+2],
            vtx_tags[i]);
  }

  fprintf(f, "Tetrahedra\n%d\n", n_tetra);
  for (int i = 0; i < n_tetra; i++) {
    fprintf(f, "%d %d %d %d %i\n",
            tetra_vtx[4*i  ],
            tetra_vtx[4*i+1],
            tetra_vtx[4*i+2],
            tetra_vtx[4*i+3],
            tetra_tag[i]);
  }

  fprintf(f, "Triangles\n%d\n", n_tri);
  for (int i = 0; i < n_tri; i++) {
    fprintf(f, "%d %d %d %i\n",
            tria_vtx[3*i    ],
            tria_vtx[3*i + 1],
            tria_vtx[3*i + 2],
            tria_tag[i]);
  }
  fprintf(f, "End\n");

  fprintf(f, "Edges\n%d\n", n_edge);
  for (int i = 0; i < n_edge; i++) {
    fprintf(f, "%d %d %i\n",
            edge_vtx[2*i    ],
            edge_vtx[2*i + 1],
            edge_tag[i]);
  }
  fprintf(f, "End\n");
  fclose(f);
}




void
PDM_write_gamma_sol
(
  const char   *filename,
  const int     n_vtx,
  const int     n_field,
  const double *fields
)
{
  PDM_UNUSED(n_field);
  // Write file
  FILE *f = fopen(filename, "w");

  fprintf(f, "MeshVersionFormatted 2\n");
  fprintf(f, "# rank %d\n\n", 0);
  fprintf(f, "Dimension\n3\n\n");
  fprintf(f, "SolAtVertices\n%d\n", n_vtx);

  fprintf(f, "%i ", n_field);
  for (int i_field = 0; i_field < n_field; i_field++) {
    fprintf(f, "1 ");
  }
  fprintf(f, "\n");
  for (int i = 0; i < n_vtx; i++) {
    for(int i_field = 0; i_field < n_field; ++i_field) {
      fprintf(f, "%20.16lf ", fields[n_field*i+i_field]);
    }
    fprintf(f, " \n");
  }
  fprintf(f, "End\n");
  fclose(f);
}


void
PDM_write_gamma_matsym
(
  const char   *filename,
  const int     n_vtx,
  const double *fields
)
{
  // Write file
  FILE *f = fopen(filename, "w");

  fprintf(f, "MeshVersionFormatted 2\n");
  fprintf(f, "# rank %d\n\n", 0);
  fprintf(f, "Dimension\n3\n\n");
  fprintf(f, "SolAtVertices\n%d\n", n_vtx);

  fprintf(f, "1 3 \n");
  for (int i = 0; i < n_vtx; i++) {
    for(int i_field = 0; i_field < 6; ++i_field) {
      fprintf(f, "%20.16lf ", fields[6*i+i_field]);
    }
    fprintf(f, " \n");
  }
  fprintf(f, "End\n");
  fclose(f);
}

/* https://pyamg.saclay.inria.fr/download/vizir/vizir4_user_guide.pdf*/
void
PDM_read_gamma_sol
(
  const char   *filename,
  const int     n_vtx,
  const int     n_field,
        double *fields
)
{
  PDM_UNUSED(n_field);

  // Read file
  FILE *f = fopen(filename, "r");

  char line[999];

  // double *lfield = fields[0];

  while (1) {

    int stat = fscanf(f, "%s", line);

    if (stat == EOF) {
      // End of file
      break;
    }

    if (strstr(line, "SolAtVertices") != NULL) {
      long _gn_vtx;
      fscanf(f, "%ld", &_gn_vtx);
      // assert(_gn_vtx == n_vtx);

      long _n_field = -1;
      long _i_kind = -1;
      fscanf(f, "%ld", &_n_field);
      for (int i_field = 0; i_field < n_field; i_field++) {
        fscanf(f, "%ld", &_i_kind);
      }

      // printf("n_field  = %i \n", (int)n_field);
      // printf("_n_field = %i \n", (int)_n_field);
      // printf("n_vtx    = %i \n", (int)n_vtx);
      // printf("_gn_vtx  = %i \n", (int)_gn_vtx);
      // printf("_i_kind  = %i \n", (int)_i_kind);

      // assert(n_field == _n_field);


      for (int i = 0; i < n_vtx; i++) {
        for (int i_field = 0; i_field < n_field; i_field++) {
          fields[n_field*i + i_field] = -(n_field*i + i_field);
          fscanf(f, "%lf", &fields[n_field*i + i_field]);
          // log_trace("fields[n_field*i + i_field] = %lf \n", fields[n_field*i + i_field]);
        }
      }
    }
  }
  fclose(f);
}
