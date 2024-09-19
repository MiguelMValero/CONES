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
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"

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
           PDM_g_num_t   *n_vtx_seg,
           double        *length,
           int           *post)
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
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
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

int main(int argc, char *argv[])
{

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                post      = 0;

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &post);

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell      = NULL;
  int          *dface_vtx_idx   = NULL;
  PDM_g_num_t  *dface_vtx       = NULL;
  double       *dvtx_coord      = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group     = NULL;
  int           dface_vtx_l;
  int           dface_group_l;

  PDM_dcube_t* dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          0.,
                                          0.,
                                          0.,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtx_l,
                         &dface_group_l);

  PDM_dcube_gen_data_get(dcube,
                          &dface_cell,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &dvtx_coord,
                          &dface_group_idx,
                          &dface_group);

  PDM_g_num_t* face_distribution = PDM_compute_entity_distribution(comm, dn_face);
  PDM_g_num_t* vtx_distribution  = PDM_compute_entity_distribution(comm, dn_vtx );
  /*
   *  Choice of extraction
   */

  PDM_g_num_t   *extract_face_distribution = NULL;
  PDM_g_num_t   *extract_vtx_distribution  = NULL;
  int           *dextract_face_vtx_idx     = NULL;
  PDM_g_num_t   *dextract_face_vtx         = NULL;
  PDM_g_num_t   *dparent_face_g_num        = NULL;
  PDM_g_num_t   *dparent_vtx_g_num         = NULL;
  PDM_g_num_t   *pextract_old_to_new       = NULL;

  PDM_dconnectivity_to_extract_dconnectivity_bis(comm,
                                                 dface_group_idx[n_face_group],
                                                 dface_group,
                                                 face_distribution,
                                                 dface_vtx_idx,
                                                 dface_vtx,
                                                 &extract_face_distribution,
                                                 &extract_vtx_distribution,
                                                 &dextract_face_vtx_idx,
                                                 &dextract_face_vtx,
                                                 &dparent_face_g_num,
                                                 &dparent_vtx_g_num,
                                                 &pextract_old_to_new);


  int dn_extract_face = extract_face_distribution[i_rank+1] - extract_face_distribution[i_rank];
  int dn_extract_vtx  = extract_vtx_distribution [i_rank+1] - extract_vtx_distribution [i_rank];

  double** tmp_dextract_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        vtx_distribution,
                                        dvtx_coord,
                                        &dn_extract_vtx,
                 (const PDM_g_num_t **) &dparent_vtx_g_num,
                                        &tmp_dextract_vtx_coord);

  double* dextract_vtx_coord = tmp_dextract_vtx_coord[0];
  free(tmp_dextract_vtx_coord);

  /*
   *  Echange de champs
   */
  int* dface_group_tag = malloc(dface_group_idx[n_face_group] * sizeof(int));

  PDM_g_num_t* dface_group_init_distrib = PDM_compute_entity_distribution(comm, dface_group_idx[n_face_group]);
  for(int i_group = 0; i_group < n_face_group; ++i_group) {
    for(int i_face = dface_group_idx[i_group]; i_face < dface_group_idx[i_group+1]; ++i_face) {
      dface_group_tag[i_face] = i_group;
    }
  }

  if(post) {
    PDM_log_trace_array_long(extract_face_distribution, n_rank+1, "extract_face_distribution:: ");
    PDM_log_trace_array_long(extract_vtx_distribution , n_rank+1, "extract_vtx_distribution::  ");

    PDM_log_trace_array_int(dextract_face_vtx_idx, dn_extract_face+1                     , "dextract_face_vtx_idx:: ");
    PDM_log_trace_array_long(dextract_face_vtx   , dextract_face_vtx_idx[dn_extract_face], "dextract_face_vtx:: "    );
    PDM_log_trace_array_long(dparent_face_g_num  , dn_extract_face                       , "dparent_face_g_num:: "   );
    PDM_log_trace_array_long(dparent_vtx_g_num   , dn_extract_vtx                        , "dparent_vtx_g_num:: "    );
  }

  /*
   *  Visulisation
   */
  PDM_g_num_t* extract_face_ln_to_gn = malloc(dn_extract_face * sizeof(PDM_g_num_t));

  for(int i = 0; i < dn_extract_face; ++i) {
    extract_face_ln_to_gn[i] = extract_face_distribution[i_rank] + i + 1;
  }

  int pn_extract_vtx = -1;
  PDM_g_num_t *pextract_vtx_ln_to_gn = NULL;
  int         *pextract_face_vtx_idx = NULL;
  int         *pextract_face_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           extract_face_distribution,
                                                           dextract_face_vtx_idx,
                                                           dextract_face_vtx,
                                                           dn_extract_face,
                                                           extract_face_ln_to_gn,
                                                           &pn_extract_vtx,
                                                           &pextract_vtx_ln_to_gn,
                                                           &pextract_face_vtx_idx,
                                                           &pextract_face_vtx);

  int** tmp_dextract_face_tag = NULL;
  PDM_part_dfield_to_pfield(comm,
                            1,
                            sizeof(int),
                            dface_group_init_distrib,
    (unsigned char    *)    dface_group_tag,
                            &dn_extract_face,
    (const PDM_g_num_t **)  &pextract_old_to_new,
    (unsigned char ***)     &tmp_dextract_face_tag);
  int* dextract_face_tag = tmp_dextract_face_tag[0];
  free(tmp_dextract_face_tag);

  /*
   * To true partition for visu
   */
  int** tmp_pextract_face_tag = NULL;
  PDM_part_dfield_to_pfield(comm,
                            1,
                            sizeof(int),
                            extract_face_distribution,
      (unsigned char    *)  dextract_face_tag,
                            &dn_extract_face,
     (const PDM_g_num_t **) &extract_face_ln_to_gn,
     (unsigned char    ***) &tmp_pextract_face_tag);
  int *pextract_face_tag = tmp_pextract_face_tag[0];
  free(tmp_pextract_face_tag);

  double** tmp_pextract_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        extract_vtx_distribution,
                                        dextract_vtx_coord,
                                        &pn_extract_vtx,
                 (const PDM_g_num_t **) &pextract_vtx_ln_to_gn,
                                        &tmp_pextract_vtx_coord);
  double* pextract_vtx_coord = tmp_pextract_vtx_coord[0];
  free(tmp_pextract_vtx_coord);

  if (post) {
    char filename[999];
    sprintf(filename, "export_face_%i.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           pn_extract_vtx,
                           pextract_vtx_coord,
                           pextract_vtx_ln_to_gn,
                           dn_extract_face,
                           pextract_face_vtx_idx,
                           pextract_face_vtx,
                           extract_face_ln_to_gn,
                           pextract_face_tag);
  }

  free(dface_group_tag);
  free(dextract_face_tag);
  free(dextract_vtx_coord);
  free(extract_face_ln_to_gn);
  free(pextract_face_tag);

  free(extract_face_distribution);
  free(extract_vtx_distribution );
  free(dextract_face_vtx_idx    );
  free(dextract_face_vtx        );
  free(dparent_face_g_num       );
  free(dparent_vtx_g_num        );
  free(pextract_old_to_new      );
  free(dface_group_init_distrib );

  free(pextract_vtx_ln_to_gn);
  free(pextract_face_vtx_idx);
  free(pextract_face_vtx    );
  free(pextract_vtx_coord   );

  free(face_distribution);
  free(vtx_distribution );

  PDM_dcube_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
